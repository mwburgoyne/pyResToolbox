"""
VLP calibration to measured well tests: fit gravity and/or friction multipliers
on a VLP method so its BHP matches measured (gauge) pressures, with an
identifiability report and leave-one-out prediction, and optionally a PI scale
on the inflow model.
"""

import math
import warnings

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from pyrestoolbox.nodal.nodal import fbhp, ipr_curve, NodalResult
from pyrestoolbox.validate import validate_choice

__all__ = ['vlp_calibrate']

_FIT_PARAMS = {'friction': ('friction_mult',), 'gravity': ('gravity_mult',),
               'both': ('gravity_mult', 'friction_mult')}
# Start grid, log-spaced over the bounds. The objective is not smooth: method
# regime switches step BHP by 0.2-8 psi, so a local search started blind can
# stall on a step; the grid picks the basin first.
_GRID_1D = 17
_GRID_2D = 9
# Leave-one-out refits start from the all-test fit instead of repeating the grid:
# dropping one test rarely moves the basin, and the grid is most of the cost. A 1-D
# refit searches this factor either side (one grid step at the default bounds) and
# falls back to the grid when its optimum lands on the window edge.
_LOO_SPAN = 1.09
_XATOL = 1e-4               # multiplier tolerance of the local refinement
_FATOL = 1e-4               # sum-of-squares tolerance (pressure units squared), 2-D refinement
# Identifiability: finite-difference step on each multiplier, relative. Large
# enough that regime-switch steps do not dominate the sensitivity.
_JAC_REL_STEP = 0.05
_CORR_WARN = 0.95           # |correlation| of the two multipliers above which they trade off
_BOUND_REL_TOL = 1e-3       # within this relative distance of a bound counts as at the bound
_RATE_KEYS = {'gas': 'qg_mscfd', 'oil': 'qt_stbpd'}


def vlp_calibrate(tests, completion, vlpmethod='WG', well_type='gas', fit='friction',
                  bounds=(0.5, 2.0), common=None, reservoir=None, ipr_kwargs=None,
                  loo=True, metric=False):
    """ Fits VLP gravity and/or friction multipliers to measured well tests.

        tests: list of dicts, one per test, each holding the measured 'bhp' (psia | barsa,
            at the bottom of the completion) and the fbhp inputs for that test, at least
            'thp' and the rate ('qg_mscfd' for gas, 'qt_stbpd' for oil), plus any that
            vary between tests (qw_bwpd, cgr, gor, wc ...)
        completion: Completion object, ending at the gauge depth
        vlpmethod: VLP method - 'HB', 'WG', 'GRAY' or 'BB'
        well_type: 'gas' or 'oil'
        fit: 'friction', 'gravity' or 'both'. The multiplier not fitted stays at 1
        bounds: (lower, upper) for each fitted multiplier, 0 < lower < upper. Default (0.5, 2.0)
        common: dict of fbhp inputs shared by every test (gsg, api, pb, rsb, gas_pvt ...)
        reservoir: Reservoir object. If given, a PI scale is also fitted: the multiplier on
            the ipr_curve rate at each measured BHP that best matches the test rates
        ipr_kwargs: dict passed to ipr_curve for the PI scale (gas_pvt, oil_pvt, wc, bo, uo, gsg)
        loo: If True (default) and there are at least two more tests than fitted
            multipliers, refit leaving each test out and predict its BHP
        metric: If True, pressures and rates in tests/common are Eclipse METRIC

        Returns NodalResult with 'gravity_mult', 'friction_mult', 'bhp_meas', 'bhp_calc',
        'residual' (calc - meas), 'rms', 'std_err' (per fitted multiplier, None with too
        few tests), 'correlation' (between the two multipliers when fit='both'),
        'at_bound' (fitted multipliers that ended on a bound), 'identifiable',
        'loo_residual' and 'loo_rms' (None when not run), 'pi_scale' and
        'ipr_rate_residual' (None without reservoir), 'n_tests'.
        Warns when the fit is not identifiable.
    """
    validate_choice(well_type, ('gas', 'oil'), 'well_type')
    validate_choice(fit, tuple(_FIT_PARAMS), 'fit')
    names = _FIT_PARAMS[fit]
    lo, hi = _check_bounds(bounds)
    tests = list(tests)
    bhp_meas, base = _fbhp_inputs(tests, dict(common or {}, completion=completion,
                                              vlpmethod=vlpmethod, well_type=well_type,
                                              metric=metric))

    def residuals(x, idx):
        mult = dict(zip(names, x))
        return np.array([fbhp(**base[i], **mult) for i in idx]) - bhp_meas[idx]

    all_idx = np.arange(len(tests))
    x = _fit(residuals, all_idx, len(names), lo, hi)
    resid = residuals(x, all_idx)

    at_bound = [n for n, v in zip(names, x)
                if v <= lo * (1 + _BOUND_REL_TOL) or v >= hi * (1 - _BOUND_REL_TOL)]
    std_err, corr = _identifiability(residuals, all_idx, x, resid)
    notes = _identifiability_notes(at_bound, bounds, corr, len(names), len(tests))
    for note in notes:
        warnings.warn(note, UserWarning, stacklevel=2)

    loo_resid = None
    if loo and len(tests) >= len(names) + 2:
        loo_resid = _leave_one_out(residuals, all_idx, x, lo, hi)

    mults = {'gravity_mult': 1.0, 'friction_mult': 1.0}
    mults.update(zip(names, (float(v) for v in x)))
    pi_scale, ipr_resid = (None, None) if reservoir is None else _pi_scale(
        tests, bhp_meas, reservoir, well_type, dict(ipr_kwargs or {}), metric)
    return NodalResult({
        **mults,
        'bhp_meas': bhp_meas, 'bhp_calc': bhp_meas + resid, 'residual': resid,
        'rms': _rms(resid),
        'std_err': std_err, 'correlation': corr, 'at_bound': at_bound,
        'identifiable': not notes,
        'loo_residual': loo_resid, 'loo_rms': _rms(loo_resid),
        'pi_scale': pi_scale, 'ipr_rate_residual': ipr_resid,
        'n_tests': len(tests),
    })


def _check_bounds(bounds):
    lo, hi = bounds
    if not (0 < lo < hi and math.isfinite(hi)):
        raise ValueError(f"bounds must satisfy 0 < lower < upper, got {bounds}")
    return lo, hi


def _fbhp_inputs(tests, shared):
    """Measured BHPs and one fbhp kwargs dict per test (shared inputs, then the test's own)."""
    if not tests:
        raise ValueError("tests must hold at least one measured test")
    bhp_meas = np.array([_measured_bhp(t, i) for i, t in enumerate(tests)])
    return bhp_meas, [dict(shared, **{k: v for k, v in t.items() if k != 'bhp'}) for t in tests]


def _rms(r):
    return None if r is None else float(np.sqrt(np.mean(r ** 2)))


def _leave_one_out(residuals, idx, x_all, lo, hi):
    """Residual of each test predicted from multipliers fitted to the others."""
    return np.array([residuals(_fit(residuals, np.delete(idx, i), len(x_all), lo, hi,
                                    x0=x_all), [i])[0]
                     for i in idx])


def _measured_bhp(test, i):
    if 'bhp' not in test or 'thp' not in test:
        raise ValueError(f"test {i} needs 'bhp' (measured) and 'thp'")
    bhp = test['bhp']
    if not (math.isfinite(bhp) and bhp > 0):
        raise ValueError(f"test {i}: measured bhp must be positive and finite, got {bhp}")
    return float(bhp)


def _identifiability_notes(at_bound, bounds, corr, n_par, n_tests):
    notes = []
    if at_bound:
        notes.append(f"{', '.join(at_bound)} ended on its bound {bounds}: the tests do not "
                     "constrain it inside the bounds")
    if corr is not None and abs(corr) > _CORR_WARN:
        notes.append(f"gravity and friction multipliers are correlated ({corr:.3f}): the tests "
                     "cannot separate them; spread the test rates or fit one multiplier")
    if n_par == 2 and n_tests < 3:
        notes.append("fitting two multipliers to fewer than 3 tests leaves no redundancy")
    return notes


def _fit(residuals, idx, n_par, lo, hi, x0=None):
    """Least-squares multipliers: log-spaced grid, then a bounded local refinement.
    With x0, the grid is skipped and the refinement starts there."""
    def sse(x):
        try:
            return float(np.sum(residuals(x, idx) ** 2))
        except RuntimeError:   # march diverged at this multiplier
            return math.inf

    if n_par == 1:
        if x0 is None:
            grid = np.geomspace(lo, hi, _GRID_1D)
            vals = [sse([g]) for g in grid]
            k = int(np.argmin(vals))
            start, f0 = grid[k], vals[k]
            a, b = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
        else:
            start, f0 = x0[0], sse(x0)
            a, b = max(lo, start / _LOO_SPAN), min(hi, start * _LOO_SPAN)
        res = minimize_scalar(lambda v: sse([v]), bounds=(a, b), method='bounded',
                              options={'xatol': _XATOL})
        edge = min(res.x - a, b - res.x) < 2 * _XATOL and lo < res.x < hi
        if x0 is not None and edge:
            return _fit(residuals, idx, n_par, lo, hi)
        return np.array([res.x if res.fun <= f0 else start])

    if x0 is None:
        grid = np.geomspace(lo, hi, _GRID_2D)
        best = min(((sse([g, f]), g, f) for g in grid for f in grid), key=lambda t: t[0])
    else:
        best = (sse(x0), x0[0], x0[1])
    res = minimize(sse, x0=[best[1], best[2]], method='Nelder-Mead', bounds=[(lo, hi)] * 2,
                   options={'xatol': _XATOL, 'fatol': _FATOL})
    return np.array(res.x if res.fun <= best[0] else best[1:])


def _identifiability(residuals, idx, x, resid):
    """Standard errors and (two-parameter) correlation from a finite-difference Jacobian."""
    n, p = len(idx), len(x)
    jac = np.empty((n, p))
    for j in range(p):
        h = _JAC_REL_STEP * x[j]
        up, dn = x.copy(), x.copy()
        up[j] += h
        dn[j] -= h
        jac[:, j] = (residuals(up, idx) - residuals(dn, idx)) / (2 * h)
    jtj = jac.T @ jac
    if np.linalg.matrix_rank(jtj) < p:
        return None, (1.0 if p == 2 else None)
    cov = np.linalg.inv(jtj)
    corr = float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1])) if p == 2 else None
    if n <= p:
        return None, corr
    s2 = float(np.sum(resid ** 2)) / (n - p)
    return [float(math.sqrt(s2 * cov[j, j])) for j in range(p)], corr


def _pi_scale(tests, bhp_meas, reservoir, well_type, ipr_kwargs, metric):
    """Least-squares multiplier s on the inflow rate: q_meas ~ s q_ipr(bhp_meas).
    Exact for rates proportional to the productivity index (Darcy, Vogel, gas
    without non-Darcy flow)."""
    key = _RATE_KEYS[well_type]
    q_meas, q_ipr = [], []
    for i, (t, bhp) in enumerate(zip(tests, bhp_meas)):
        if key not in t:
            raise ValueError(f"test {i} needs '{key}' for the PI scale")
        kw = dict(ipr_kwargs)
        if well_type == 'oil' and 'wc' in t:
            kw['wc'] = t['wc']
        q_meas.append(t[key])
        q_ipr.append(ipr_curve(reservoir, well_type=well_type, n_points=2, min_pwf=bhp,
                               metric=metric, **kw)['rate'][0])
    q_meas, q_ipr = np.array(q_meas, dtype=float), np.array(q_ipr, dtype=float)
    denom = float(np.sum(q_ipr ** 2))
    if denom <= 0:
        raise ValueError("the inflow model gives zero rate at every measured BHP")
    scale = float(np.sum(q_ipr * q_meas)) / denom
    return scale, scale * q_ipr - q_meas
