"""
Capacitance-resistance models (CRM) for waterfloods: CRMT (field tank) and CRMP
(one time constant per producer), with optional producer BHP, after

Sayarpour, Zuluaga, Kabir & Lake (2009), The use of capacitance-resistance models for
rapid estimation of waterflood performance and optimization, JPSE 69, 227-238
(doi 10.1016/j.petrol.2009.09.006). Equation numbers below are that paper's.

Producer j obeys dq_j/dt + q_j / tau_j = (1 / tau_j) sum_i f_ij i_i(t) - J_j dpwf_j/dt
(Eq. 1). With injection constant and BHP linear over each interval the exact solution
(Eq. 5) is the one-step recursion

    q_j(t_k) = q_j(t_{k-1}) a + (1 - a) [sum_i f_ij I_i(k) - J_j tau_j dpwf_j(k) / dt_k],
    a = exp(-dt_k / tau_j).

Unit-agnostic: rates in any consistent reservoir-volume rate unit (injection and
total-liquid production in the same unit), time in any unit (tau comes back in it), BHP
in any pressure unit (J in rate per pressure).
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

from pyrestoolbox.shared_fns import validate_pe_inputs
from pyrestoolbox.validate import validate_choice

__all__ = ['crm_fit', 'crm_oil_cut_fit', 'CRMResult']

# Bounds on the fitted time constant, as multiples of the shortest step and of the
# history span
_TAU_MIN_FACTOR = 0.01
_TAU_MAX_FACTOR = 100.0
# Log-spaced common-tau grid used to initialise the joint fit
_TAU_GRID_POINTS = 25
# Slack on sum_j f_ij <= 1 before the constrained polish runs
_SUM_F_TOL = 1e-9


@dataclass
class CRMResult:
    """Fitted capacitance-resistance model.

    model: 'CRMP' or 'CRMT'
    f: Connectivities, shape (n_injectors, n_producers); f[i, j] is the steady-state
        fraction of injector i's rate produced by producer j (CRMT: [[1.0]])
    tau: Time constants per producer (time units of t)
    q0: Fitted rate at t[0] per producer
    J: Productivity index per producer (rate per pressure), None without BHP
    r_squared: Per-producer coefficient of determination of the liquid-rate match
    q_fit: Fitted rates, shape (n_times, n_producers)
    t: History times
    """
    model: str
    f: np.ndarray
    tau: np.ndarray
    q0: np.ndarray
    J: Optional[np.ndarray]
    r_squared: np.ndarray
    q_fit: np.ndarray
    t: np.ndarray

    def predict(self, t, inj, bhp=None):
        """Producer rates over a schedule starting at the history's t[0] (pass the
        history plus any future period). inj shaped like the fitting input; bhp
        required when the model was fitted with BHP. Returns (n_times, n_producers),
        or a 1-D array for CRMT."""
        t, inj, _, bhp = _shape_inputs(t, inj, None, bhp, self.model)
        if (self.J is None) != (bhp is None):
            raise ValueError("bhp must be given exactly when the model was fitted with it")
        if t[0] != self.t[0]:
            raise ValueError("t must start at the fitted history's first time")
        q = _simulate(np.diff(t), inj, _bhp_steps(bhp), self.f, self.tau, self.q0, self.J)
        return q[:, 0] if self.model == 'CRMT' else q

    def pore_volume(self, ct):
        """Drainage pore volume per producer from tau = ct Vp / J (Eq. 2), in rate x
        time units (rb/d x days = rb). Needs J, so only for fits with BHP."""
        if self.J is None:
            raise ValueError("Pore volume needs J; fit with bhp to obtain it")
        if ct <= 0:
            raise ValueError(f"ct must be positive, got {ct}")
        return self.tau * self.J / ct


def _bhp_steps(bhp):
    """BHP change over each interval, shape (n-1, n_producers), or None."""
    return None if bhp is None else np.diff(bhp, axis=0)


def _simulate(dt, inj, dbhp, f, tau, q0, J):
    """Eq. (5) recursion. inj[k] is the injection over (t[k-1], t[k]]; inj[0] unused."""
    a = np.exp(-dt[:, None] / tau[None, :])
    drive = inj[1:] @ f
    if J is not None:
        drive = drive - J * tau * dbhp / dt[:, None]
    q = np.empty((len(dt) + 1, len(tau)))
    q[0] = q0
    for k in range(len(dt)):
        q[k + 1] = q[k] * a[k] + (1.0 - a[k]) * drive[k]
    return q


def _shape_inputs(t, inj, prod, bhp, model):
    """Arrays as (n,), (n, Ni), (n, Np), (n, Np) or None; CRMT lumps to one well."""
    t = np.asarray(t, dtype=float)
    inj = np.asarray(inj, dtype=float)
    inj = inj[:, None] if inj.ndim == 1 else inj
    if prod is not None:
        prod = np.asarray(prod, dtype=float)
        prod = prod[:, None] if prod.ndim == 1 else prod
    if bhp is not None:
        bhp = np.asarray(bhp, dtype=float)
        bhp = bhp[:, None] if bhp.ndim == 1 else bhp
    if model == 'CRMT':
        inj = inj.sum(axis=1, keepdims=True)
        prod = None if prod is None else prod.sum(axis=1, keepdims=True)
    return t, inj, prod, bhp


def _check_inputs(t, inj, prod, bhp):
    n = len(t)
    if t.ndim != 1 or n < 3 or np.any(np.diff(t) <= 0):
        raise ValueError("t must be 1-D, strictly increasing, with at least 3 points")
    if inj.shape[0] != n or prod.shape[0] != n:
        raise ValueError(f"inj and prod need {n} rows (one per time), got {inj.shape[0]} "
                         f"and {prod.shape[0]}")
    for name, arr in (('inj', inj), ('prod', prod)):
        if not np.all(np.isfinite(arr)) or np.any(arr < 0):
            raise ValueError(f"{name} must be finite and non-negative")
    if bhp is not None:
        if bhp.shape != prod.shape:
            raise ValueError(f"bhp must have the shape of prod {prod.shape}, got {bhp.shape}")
        validate_pe_inputs(p=bhp)


def _filter(dt, u, tau):
    """Zero-initial-rate response to the drives u, shape (n-1, m), for one tau:
    y_k = a y_{k-1} + (1 - a) u_k, the Eq. (5) recursion."""
    a = np.exp(-dt / tau)[:, None]
    y = np.zeros((len(dt) + 1, u.shape[1]))
    for k in range(len(dt)):
        y[k + 1] = y[k] * a[k] + (1.0 - a[k]) * u[k]
    return y


class _Producer:
    """One producer's least squares at fixed tau. The model is linear in
    (f_1j..f_Nij, q0_j, J_j), so for each tau bounded least squares gives the exact
    optimum; the remaining search is one-dimensional in tau (variable projection)."""

    def __init__(self, t, inj, y, dbhp_j, fit_f):
        self.t, self.dt, self.inj, self.y = t, np.diff(t), inj, y
        self.dbhp, self.fit_f = dbhp_j, fit_f
        self.ni = inj.shape[1] if fit_f else 0

    def design(self, tau):
        """Columns and the fixed part of the target (CRMT: f = 1 is not fitted)."""
        x = _filter(self.dt, self.inj[1:], tau)
        cols = [x] if self.fit_f else []
        cols.append(np.exp(-(self.t - self.t[0]) / tau)[:, None])
        if self.dbhp is not None:
            cols.append(-_filter(self.dt, (tau * self.dbhp / self.dt)[:, None], tau))
        target = self.y if self.fit_f else self.y - x[:, 0]
        return np.hstack(cols), target

    def solve(self, tau, lam):
        """theta minimising |A theta - target|^2 + lam . f under the bounds, and that
        objective. lam (n_inj,) are the Eq. (20) multipliers, None when not active."""
        from scipy.optimize import lsq_linear
        a, target = self.design(tau)
        upper = np.r_[np.ones(self.ni), np.full(a.shape[1] - self.ni, np.inf)]
        shifted = target
        if lam is not None and np.any(lam > 0):
            c = np.r_[lam, np.zeros(a.shape[1] - self.ni)]
            shifted = target - a @ np.linalg.lstsq(a.T @ a, 0.5 * c, rcond=None)[0]
        theta = lsq_linear(a, shifted, bounds=(0.0, upper), method='bvls').x
        obj = float(np.sum((a @ theta - target) ** 2))
        if lam is not None:
            obj += float(lam @ theta[:self.ni])
        return theta, obj

    def fit(self, taus, lam=None, hint=None):
        """Best tau: the grid minimum refined by bounded Brent search in ln tau. With
        a hint (the previous tau, inside the dual iterations) the grid is skipped and
        Brent searches one grid step either side of it."""
        from scipy.optimize import minimize_scalar
        step = np.log(taus[1] / taus[0])
        if hint is None:
            objs = [self.solve(tau, lam)[1] for tau in taus]
            i = int(np.argmin(objs))
            best_tau, best_obj = float(taus[i]), objs[i]
        else:
            best_tau, best_obj = hint, self.solve(hint, lam)[1]
        lo = max(np.log(best_tau) - step, np.log(taus[0]))
        hi = min(np.log(best_tau) + step, np.log(taus[-1]))
        res = minimize_scalar(lambda v: self.solve(np.exp(v), lam)[1], bounds=(lo, hi),
                              method='bounded', options={'xatol': 1e-7})
        tau = float(np.exp(res.x)) if res.fun < best_obj else best_tau
        theta, obj = self.solve(tau, lam)
        return tau, theta, obj


def _fit_all(producers, taus, lam=None, hints=None):
    hints = hints or [None] * len(producers)
    return [p.fit(taus, lam, h) for p, h in zip(producers, hints)]


def _coupled_fit(producers, taus, ni, scale, fits):
    """Enforce sum_j f_ij <= 1 (Eq. 20) through one Lagrange multiplier per injector.
    For given multipliers the Lagrangian separates by producer; the dual is maximised
    over lam >= 0 with gradient sum_j f_ij - 1 (Danskin), lam scaled by `scale`."""
    from scipy.optimize import minimize

    hints = [fit[0] for fit in fits]

    def neg_dual(mu):
        trial = _fit_all(producers, taus, mu * scale, hints)
        hints[:] = [fit[0] for fit in trial]
        f_sum = np.sum([fit[1][:ni] for fit in trial], axis=0)
        dual = sum(fit[2] for fit in trial) - float(mu * scale @ np.ones(ni))
        return -dual / scale, -(f_sum - 1.0)

    mu = minimize(neg_dual, np.zeros(ni), jac=True, method='L-BFGS-B',
                  bounds=[(0.0, None)] * ni, options={'gtol': 1e-10}).x
    # Final pass from the full grid, so the warm starts cannot have trapped a producer
    return _fit_all(producers, taus, mu * scale)


def _feasible_f(fits, producers, ni):
    """Scale any injector's f row whose sum still exceeds one (duality gap or
    tolerance) down to one, then refit q0 and J with tau and f fixed."""
    from scipy.optimize import lsq_linear
    f = np.array([fit[1][:ni] for fit in fits]).T
    f = f / np.maximum(f.sum(axis=1), 1.0)[:, None]
    out = []
    for j, (tau, theta, _) in enumerate(fits):
        a, target = producers[j].design(tau)
        rest = lsq_linear(a[:, ni:], target - a[:, :ni] @ f[:, j], bounds=(0.0, np.inf),
                          method='bvls').x
        out.append((tau, np.r_[f[:, j], rest]))
    return out


def crm_fit(t, inj, prod, bhp=None, model='CRMP'):
    """Fit a capacitance-resistance model to injection and total-liquid production.

    Least squares with f_ij in [0, 1], q0 >= 0, J >= 0 and, for each injector, the
    sum over producers of f_ij <= 1 (Eq. 20). For fixed tau the model is linear in
    the other parameters, so each producer is solved exactly by bounded least squares
    inside a one-dimensional tau search (grid, then Brent in ln tau). When Eq. (20)
    binds, one Lagrange multiplier per injector couples the producers and the dual is
    maximised; any residual excess is scaled out. Sayarpour et al. minimise absolute
    error and
    recommend at least four data points per unknown: CRMP has n_inj + 2 unknowns per
    producer (+1 with BHP), CRMT two. CRMT lumps all wells into one tank with f = 1
    (Eq. 8); BHP is not used for CRMT.

    t: Times, strictly increasing (any unit; tau is returned in it)
    inj: Injection rates, shape (n_times,) or (n_times, n_injectors), reservoir
        volume rate; inj[k] is the average over (t[k-1], t[k]] (inj[0] is unused)
    prod: Total-liquid production rates, shape (n_times,) or (n_times, n_producers),
        same rate unit as inj, at each time
    bhp: Optional producer bottomhole pressures, same shape as prod (any pressure
        unit), adding a productivity index J per producer to the fit
    model: 'CRMP' (default) or 'CRMT'

    Returns CRMResult with f, tau, q0, J, per-producer r_squared, the fitted rates,
    and predict() / pore_volume() methods.
    """
    validate_choice(model, ('CRMP', 'CRMT'), 'model')
    if model == 'CRMT':
        bhp = None
    t, inj, prod, bhp = _shape_inputs(t, inj, prod, bhp, model)
    _check_inputs(t, inj, prod, bhp)

    dbhp = _bhp_steps(bhp)
    dt, span = np.diff(t), t[-1] - t[0]
    tau_bounds = (_TAU_MIN_FACTOR * dt.min(), _TAU_MAX_FACTOR * span)
    taus = np.geomspace(*tau_bounds, _TAU_GRID_POINTS)
    fit_f = model == 'CRMP'
    ni = inj.shape[1] if fit_f else 0
    producers = [_Producer(t, inj, prod[:, j], None if dbhp is None else dbhp[:, j], fit_f)
                 for j in range(prod.shape[1])]
    fits = _fit_all(producers, taus)
    if fit_f and np.sum([fit[1][:ni] for fit in fits], axis=0).max() > 1.0 + _SUM_F_TOL:
        fits = _coupled_fit(producers, taus, ni, float(np.mean(prod)) ** 2 * len(t), fits)
        fits = _feasible_f(fits, producers, ni)
    tau = np.array([fit[0] for fit in fits])
    theta = np.array([fit[1] for fit in fits])            # (n_producers, n_params)
    f = theta[:, :ni].T if fit_f else np.ones((1, 1))
    q0 = theta[:, ni]
    J = theta[:, ni + 1] if bhp is not None else None
    q_fit = _simulate(dt, inj, dbhp, f, tau, q0, J)
    ss_res = np.sum((prod - q_fit) ** 2, axis=0)
    ss_tot = np.sum((prod - prod.mean(axis=0)) ** 2, axis=0)
    r2 = np.where(ss_tot > 0, 1.0 - ss_res / np.where(ss_tot > 0, ss_tot, 1.0), np.nan)
    return CRMResult(model=model, f=f, tau=tau, q0=q0, J=J, r_squared=r2, q_fit=q_fit, t=t)


def crm_oil_cut_fit(cwi, fo):
    """Fit the power-law oil-cut model f_o = 1 / (1 + alpha CWI^beta) (Eq. 15).

    ln(WOR) = ln(alpha) + beta ln(CWI) is fitted by least squares, using the points
    with 0 < fo < 1 and cwi > 0. Sayarpour et al. caution the model holds only on the
    straight-line part of a log-log WOR vs cumulative-injection plot (mature,
    high-watercut floods), so pass only that section. For a CRMP producer, CWI is
    its connected cumulative injection, the time integral of sum_i f_ij i_i (Eq. 18).

    cwi: Cumulative water injected (any volume unit; alpha depends on it)
    fo: Oil cut (fraction of total liquid) at each cwi

    Returns dict with 'alpha', 'beta', 'r_squared' (in log WOR) and 'n_points'.
    """
    cwi = np.asarray(cwi, dtype=float)
    fo = np.asarray(fo, dtype=float)
    if cwi.shape != fo.shape or cwi.ndim != 1:
        raise ValueError("cwi and fo must be 1-D arrays of the same length")
    keep = (cwi > 0) & (fo > 0) & (fo < 1) & np.isfinite(cwi) & np.isfinite(fo)
    if keep.sum() < 3:
        raise ValueError("Need at least 3 points with cwi > 0 and 0 < fo < 1")
    x, y = np.log(cwi[keep]), np.log(1.0 / fo[keep] - 1.0)
    beta, ln_alpha = np.polyfit(x, y, 1)
    resid = y - (ln_alpha + beta * x)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan
    return {'alpha': float(np.exp(ln_alpha)), 'beta': float(beta),
            'r_squared': float(r2), 'n_points': int(keep.sum())}
