#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
    pyResToolbox - A collection of Reservoir Engineering Utilities
              Copyright (C) 2022, Mark Burgoyne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    The GNU General Public License can be found in the LICENSE directory,
    and at  <https://www.gnu.org/licenses/>.

          Contact author at mark.w.burgoyne@gmail.com

Material balance calculations for gas and oil reservoirs.

Functions
---------
gas_matbal      P/Z gas material balance (OGIP estimation)
oil_matbal      Havlena-Odeh oil material balance (OOIP estimation)

Classes
-------
GasMatbalResult     Result from gas_matbal()
OilMatbalResult     Result from oil_matbal()
"""

__all__ = [
    'gas_matbal', 'oil_matbal',
    'GasMatbalResult', 'OilMatbalResult',
]

import warnings

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List
from scipy.optimize import least_squares

import pyrestoolbox.gas as gas
import pyrestoolbox.oil as oil
from pyrestoolbox.constants import CUFTperBBL, psc, tscr, degF2R
from pyrestoolbox.constants import BAR_TO_PSI, degc_to_degf, PSI_TO_BAR
from pyrestoolbox.constants import SCF_PER_STB_TO_SM3_PER_SM3, SM3_PER_SM3_TO_SCF_PER_STB
from pyrestoolbox.classes import z_method, c_method, rs_method, bo_method
from pyrestoolbox.shared_fns import ransac_linreg


@dataclass
class GasMatbalResult:
    """Result from gas material balance.

    Attributes
    ----------
    ogip : float
        Original gas in place (same volume units as Gp input).
    pz : np.ndarray
        P/Z values at each pressure step.
    gp : np.ndarray
        Cumulative gas production at each step.
    slope : float
        Slope of the P/Z vs Gp regression line.
    intercept : float
        Intercept of the P/Z vs Gp regression line.
    r_squared : float
        Coefficient of determination of the linear fit.
    p_initial : float
        Initial reservoir pressure.
    z_initial : float
        Z-factor at initial pressure.
    bg : np.ndarray
        Gas FVF (rcf/scf) at each pressure step.
    F : np.ndarray
        Underground withdrawal at each step.
    Et : np.ndarray
        Gas expansion term (Bg - Bgi) at each step.
    cole_F_over_Et : np.ndarray
        Cole plot diagnostic: F/Et at each step (NaN at index 0).
    method : str
        OGIP determination method: 'pz' or 'havlena_odeh'.
    pz_curvature : float
        Signed bend of P/Z vs Gp: the midpoint offset of a quadratic fit from its
        chord, as a fraction of the P/Z decline over the data. Positive (concave
        up) is consistent with pressure support (aquifer influx, a second
        compartment); negative with abnormal pressure (compaction) or a growing
        drainage area. Either way a straight-line OGIP is unreliable when it is
        large compared with the scatter. NaN with fewer than four points.
    """
    ogip: float
    pz: np.ndarray
    gp: np.ndarray
    slope: float
    intercept: float
    r_squared: float
    p_initial: float
    z_initial: float
    bg: np.ndarray = field(default_factory=lambda: np.array([]))
    F: np.ndarray = field(default_factory=lambda: np.array([]))
    Et: np.ndarray = field(default_factory=lambda: np.array([]))
    cole_F_over_Et: np.ndarray = field(default_factory=lambda: np.array([]))
    method: str = 'pz'
    pz_curvature: float = float('nan')


def _pz_curvature(gp, pz):
    """Signed midpoint offset of a quadratic P/Z fit from its chord, per unit P/Z decline."""
    if len(gp) < 4 or np.ptp(gp) <= 0:
        return float('nan')
    c2, c1, c0 = np.polyfit(gp, pz, 2)
    g0, g1 = gp.min(), gp.max()
    fall = abs((c2 * g1 ** 2 + c1 * g1) - (c2 * g0 ** 2 + c1 * g0))
    if fall <= 0:
        return float('nan')
    return float(c2 * (g1 - g0) ** 2 / 4.0 / fall)


@dataclass
class OilMatbalResult:
    """Result from oil material balance (Havlena-Odeh).

    Attributes
    ----------
    ooip : float
        Original oil in place (STB | sm3).
    F : np.ndarray
        Underground withdrawal at each step.
    Eo : np.ndarray
        Oil expansion term at each step.
    Eg : np.ndarray
        Gas cap expansion term at each step.
    Efw : np.ndarray
        Formation and water compressibility term at each step.
    drive_indices : dict
        Drive index fractions at each step from the fitted OOIP: 'DDI' (depletion),
        'SDI' (segregation/gas cap), 'CDI' (compaction/expansion), 'WDI' (water/aquifer
        influx), each N x term / F (np.ndarray, 0 where F is 0). They sum to 1 only where
        the material balance closes, so the sum is a diagnostic. 'WDI' is all-zero when no
        ``We`` aquifer influx is supplied.
    p : np.ndarray
        Pressure array used.
    pvt : dict
        PVT properties at each step: 'Rs', 'Bo', 'Bg' arrays.
    regressed : dict or None
        Regressed parameter values, in the caller's units, when ``regress`` is used.
    ooip_points : np.ndarray
        Per-point OOIP (F - We)/Et (NaN where Et is 0). Plotted against F it is the
        Campbell plot: a flat trend supports the model, a rising one unaccounted support.
    regress_report : dict or None
        With ``regress``: 'at_bound' (parameters that ended on a bound), 'max_correlation'
        (largest |correlation| between OOIP and the regressed parameters), 'identifiable'.
    """
    ooip: float
    F: np.ndarray
    Eo: np.ndarray
    Eg: np.ndarray
    Efw: np.ndarray
    drive_indices: dict
    p: np.ndarray
    pvt: dict
    regressed: Optional[dict] = None
    ooip_points: np.ndarray = field(default_factory=lambda: np.array([]))
    regress_report: Optional[dict] = None


def gas_matbal(p, Gp, degf, sg=0.65, co2=0, h2s=0, n2=0, h2=0,
               Wp=None, Bw=1.0, We=None,
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None, he=0, ag=False):
    """P/Z gas material balance for OGIP estimation.

    Performs linear regression of P/Z vs cumulative gas production to
    determine original gas in place (OGIP = -intercept/slope). Optionally
    computes Cole plot diagnostics (F/Et vs Gp) and Havlena-Odeh
    regression when cumulative water influx (We) is provided.

    Parameters
    ----------
    p : array-like
        Reservoir pressures at each survey (psia | barsa). First value is
        initial pressure (Gp=0 at that point, or the first Gp entry).
    Gp : array-like
        Cumulative gas production at each pressure survey. Same length as p.
        Units are user-defined (e.g. Bscf, MMscf) — OGIP will be in the same units.
        When Wp or We are provided, Gp should be in scf (or sm3 if metric)
        for dimensional consistency with Bg.
    degf : float
        Reservoir temperature (deg F | deg C).
    sg : float
        Gas specific gravity (default 0.65).
    co2 : float
        CO2 mole fraction (default 0).
    h2s : float
        H2S mole fraction (default 0).
    n2 : float
        N2 mole fraction (default 0).
    h2 : float
        H2 mole fraction (default 0).
    he : float
        He mole fraction (default 0; keyword, last in the signature). Like H2,
        it forces the BNS Z-factor.
    Wp : array-like, optional
        Cumulative water production (STB | sm3). Same length as p.
    Bw : float
        Water FVF (rb/stb | rm3/sm3, default 1.0). Used when Wp provided.
    We : array-like, optional
        Cumulative water influx (rcf | rm3). Same length as p. When
        provided, Havlena-Odeh regression is used for OGIP instead of P/Z.
        aquifer_influx returns rb in field units: multiply by 5.61458 for rcf.
    zmethod : str or z_method
        Z-factor method (default 'DAK').
    cmethod : str or c_method
        Critical property method (default 'PMC').
    metric : bool
        If True, p in barsa and degf in deg C (default False).
    pvt_table : dict, optional
        Tabulated gas PVT. Either ``{'p': [...], 'Z': [...]}`` or
        ``{'p': [...], 'Bg': [...]}``. Pressures in psia|barsa, Bg in
        rcf/scf|rm3/sm3. Providing both 'Z' and 'Bg' keys raises ValueError.

    Returns
    -------
    GasMatbalResult
    """
    p = np.asarray(p, dtype=float)
    Gp = np.asarray(Gp, dtype=float)

    if np.any(p <= 0):
        raise ValueError("Pressures must be positive")
    if len(p) != len(Gp):
        raise ValueError(f"p and Gp must have same length, got {len(p)} and {len(Gp)}")
    if len(p) < 2:
        raise ValueError("Need at least 2 pressure/production data points")

    n = len(p)

    # Validate optional arrays
    if Wp is not None:
        Wp = np.asarray(Wp, dtype=float)
        if len(Wp) != n:
            raise ValueError(f"Wp must have same length as p, got {len(Wp)} and {n}")
    if We is not None:
        We = np.asarray(We, dtype=float)
        if len(We) != n:
            raise ValueError(f"We must have same length as p, got {len(We)} and {n}")

    if pvt_table is not None:
        # Validate pvt_table
        has_z = 'Z' in pvt_table
        has_bg = 'Bg' in pvt_table
        if has_z and has_bg:
            raise ValueError("pvt_table must contain either 'Z' or 'Bg', not both")
        if not has_z and not has_bg:
            raise ValueError("pvt_table must contain 'p' and either 'Z' or 'Bg'")
        if 'p' not in pvt_table:
            raise ValueError("pvt_table must contain 'p' key")

        pvt_p = np.asarray(pvt_table['p'], dtype=float)
        sort_idx = np.argsort(pvt_p)
        pvt_p = pvt_p[sort_idx]

        # Convert metric pressures for internal use
        if metric:
            degf_field = degc_to_degf(degf)
            p_psia = p * BAR_TO_PSI
            pvt_p_psia = pvt_p * BAR_TO_PSI
        else:
            degf_field = degf
            p_psia = p
            pvt_p_psia = pvt_p

        degR = degf_field + degF2R

        # Check survey pressures within table range
        if np.any(p < pvt_p[0] - 1e-6) or np.any(p > pvt_p[-1] + 1e-6):
            raise ValueError("Survey pressures must fall within pvt_table pressure range")

        if has_z:
            pvt_Z = np.asarray(pvt_table['Z'], dtype=float)[sort_idx]
            z = np.interp(p, pvt_p, pvt_Z)
            # Compute Bg from Z: Bg = Z * T_R / (p_psia * tscr/psc)  (rcf/scf)
            bg = z * degR / (p_psia * (tscr / psc))
        else:
            pvt_Bg = np.asarray(pvt_table['Bg'], dtype=float)[sort_idx]
            # rcf/scf and rm3/sm3 are numerically identical — no conversion
            bg = np.interp(p, pvt_p, pvt_Bg)
            # Back-compute Z from Bg
            z = bg * p_psia * (tscr / psc) / degR
    else:
        # Z-factor at each pressure (gas_z handles metric internally)
        z = np.array([
            gas.gas_z(pi, sg, degf, co2=co2, h2s=h2s, n2=n2, h2=h2, he=he, ag=ag,
                      zmethod=zmethod, cmethod=cmethod, metric=metric)
            for pi in p
        ])

        # Bg at each pressure (rcf/scf — gas_bg handles metric internally)
        bg = np.array([
            gas.gas_bg(pi, sg, degf, co2=co2, h2s=h2s, n2=n2, h2=h2, he=he, ag=ag,
                       zmethod=zmethod, cmethod=cmethod, metric=metric)
            for pi in p
        ])

    pz = p / z

    # Linear regression: P/Z = intercept + slope * Gp (RANSAC for outlier robustness)
    slope, intercept, _ = ransac_linreg(Gp, pz)

    # P/Z OGIP = -intercept / slope (where P/Z = 0)
    if abs(slope) < 1e-30:
        raise RuntimeError("P/Z vs Gp regression has near-zero slope — cannot determine OGIP")
    pz_ogip = -intercept / slope

    # R-squared
    pz_pred = slope * Gp + intercept
    ss_res = np.sum((pz - pz_pred) ** 2)
    ss_tot = np.sum((pz - np.mean(pz)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Gas expansion: Et = Bg - Bgi
    Et = bg - bg[0]

    # Underground withdrawal: F = Gp * Bg (+ water term if Wp provided)
    F = Gp * bg
    if Wp is not None:
        water_factor = 1.0 if metric else CUFTperBBL
        F = F + Wp * Bw * water_factor

    # Cole plot diagnostic: F / Et (NaN at index 0 where Et = 0)
    cole_F_over_Et = np.full(n, np.nan)
    nonzero_Et = np.abs(Et) > 1e-30
    cole_F_over_Et[nonzero_Et] = F[nonzero_Et] / Et[nonzero_Et]

    # OGIP determination
    if We is not None:
        # Havlena-Odeh: forced-through-origin regression
        # (F - We) = G * Et  =>  G = sum((F-We)*Et) / sum(Et^2)
        valid = np.abs(Et) > 1e-30
        if not np.any(valid):
            raise RuntimeError("Cannot compute Havlena-Odeh OGIP — all Et values are zero")
        F_net = F[valid] - We[valid]
        Et_valid = Et[valid]
        ogip = float(np.sum(F_net * Et_valid) / np.sum(Et_valid ** 2))
        method = 'havlena_odeh'
    else:
        ogip = pz_ogip
        method = 'pz'

    return GasMatbalResult(
        ogip=float(ogip),
        pz=pz,
        gp=Gp,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
        p_initial=float(p[0]),
        z_initial=float(z[0]),
        bg=bg,
        F=F,
        Et=Et,
        cole_F_over_Et=cole_F_over_Et,
        method=method,
        pz_curvature=_pz_curvature(Gp, pz),
    )


_ALLOWED_REGRESS = {'m', 'cf', 'cw', 'sw_i'}
# Regression identifiability: a parameter within this fraction of its bound range of
# a bound is reported at the bound; |correlation| between OOIP and a parameter (or two
# parameters) above _MB_CORR_WARN means the data fix only a combination of them.
_MB_BOUND_TOL = 1e-3
_MB_CORR_WARN = 0.99
_MB_JAC_STEP = 1e-3         # finite-difference step, fraction of each bound range


def _havlena_odeh_et(Eo, Eg, dp, Boi, m, cf, cw, sw_i):
    """Total expansion Et = Eo + m Eg + (1 + m) Efw (Havlena and Odeh 1963) and Efw."""
    Efw = Boi * (cw * sw_i + cf) / (1.0 - sw_i) * dp if sw_i < 1.0 else np.zeros_like(dp)
    return Eo + m * Eg + (1.0 + m) * Efw, Efw


def _ooip_slope(F_net, Et):
    """Least-squares slope of F - We = N Et through the origin."""
    ss = float(np.dot(Et, Et))
    return float(np.dot(Et, F_net)) / ss if ss > 0 else None


def _regress_oil(regress, base_vals, F_net, Eo, Eg, dp, Boi, metric):
    """Least squares on F - We - N Et over the regressed parameters, N profiled out
    (the slope for each trial set). Returns the fitted values (internal units) and a
    report of parameters on a bound and correlations from the Jacobian."""
    names = list(regress)
    lo = np.array([regress[k][0] for k in names], dtype=float)
    hi = np.array([regress[k][1] for k in names], dtype=float)
    if metric:
        conv = np.array([1.0 / BAR_TO_PSI if k in ('cw', 'cf') else 1.0 for k in names])
        lo, hi = lo * conv, hi * conv
    span = hi - lo
    scale = max(float(np.max(np.abs(F_net))), 1e-30)

    def model(u):
        vals = dict(base_vals, **dict(zip(names, lo + u * span)))
        Et, _ = _havlena_odeh_et(Eo, Eg, dp, Boi, vals['m'], vals['cf'], vals['cw'], vals['sw_i'])
        return Et

    def resid(u):
        Et = model(u)
        N = _ooip_slope(F_net, Et)
        return (F_net - (N or 0.0) * Et) / scale

    x0 = np.array([base_vals[k] for k in names], dtype=float)
    u0 = np.where((x0 >= lo) & (x0 <= hi), (x0 - lo) / span, 0.5)
    u = least_squares(resid, u0, bounds=(0.0, 1.0)).x
    fitted = dict(zip(names, lo + u * span))
    return fitted, _regress_report(names, u, model, F_net)


def _regress_report(names, u, model, F_net):
    """Parameters on a bound, and the largest |correlation| among N and the regressed
    parameters from the column-normalised finite-difference Jacobian of N Et."""
    at_bound = [k for k, v in zip(names, u) if min(v, 1.0 - v) < _MB_BOUND_TOL]
    Et = model(u)
    N = _ooip_slope(F_net, Et) or 0.0
    cols = [Et]
    for j in range(len(u)):
        up, dn = u.copy(), u.copy()
        up[j] += _MB_JAC_STEP
        dn[j] -= _MB_JAC_STEP
        cols.append(N * (model(up) - model(dn)) / (2 * _MB_JAC_STEP))
    jac = np.column_stack(cols)
    norms = np.linalg.norm(jac, axis=0)
    if np.any(norms <= 0) or np.linalg.matrix_rank(jac / np.where(norms > 0, norms, 1.0)) < jac.shape[1]:
        max_corr = 1.0
    else:
        cov = np.linalg.inv((jac / norms).T @ (jac / norms))
        with np.errstate(invalid='ignore'):
            d = np.sqrt(np.diag(cov))
        corr = cov / np.outer(d, d)
        max_corr = float(np.max(np.abs(corr - np.eye(len(d)))))
        if not np.isfinite(max_corr):   # numerically collinear: inverse lost definiteness
            max_corr = 1.0
    notes = []
    if at_bound:
        notes.append(f"regressed {', '.join(at_bound)} ended on its bound: the data do not "
                     "constrain it inside the bounds")
    if max_corr > _MB_CORR_WARN:
        notes.append(f"OOIP and the regressed parameters are correlated ({max_corr:.4f}): the "
                     "data fix only a combination of them (an undersaturated history cannot "
                     "separate N from cf); fix one or add data below the bubble point")
    for note in notes:
        warnings.warn(note, UserWarning, stacklevel=3)
    return {'at_bound': at_bound, 'max_correlation': max_corr, 'identifiable': not notes}


def _validate_regress(regress):
    """Check the oil_matbal regress dict: known keys, (lower, upper) bounds."""
    if regress is None:
        return
    for key, val in regress.items():
        if key not in _ALLOWED_REGRESS:
            raise ValueError(f"Invalid regress key '{key}'. Allowed: {sorted(_ALLOWED_REGRESS)}")
        if not isinstance(val, (tuple, list)) or len(val) != 2:
            raise ValueError(f"regress['{key}'] must be a (lower, upper) tuple, got {val}")
        if val[0] >= val[1]:
            raise ValueError(f"regress['{key}'] bounds must have lower < upper, got {val}")


def _prepare_pvt_table(pvt_table, metric):
    """Validate a tabulated oil PVT dict and return (p, Rs, Bo, Bg) sorted by
    pressure in oilfield units (psia, scf/stb, rb/stb, rb/scf)."""
    required = {'p', 'Rs', 'Bo', 'Bg'}
    missing = required - set(pvt_table.keys())
    if missing:
        raise ValueError(f"pvt_table missing required keys: {missing}")
    pvt_p = np.asarray(pvt_table['p'], dtype=float)
    pvt_Rs = np.asarray(pvt_table['Rs'], dtype=float)
    pvt_Bo = np.asarray(pvt_table['Bo'], dtype=float)
    pvt_Bg = np.asarray(pvt_table['Bg'], dtype=float)
    if not (len(pvt_p) == len(pvt_Rs) == len(pvt_Bo) == len(pvt_Bg)):
        raise ValueError("All pvt_table arrays must have the same length")

    sort_idx = np.argsort(pvt_p)
    pvt_p = pvt_p[sort_idx]
    pvt_Rs = pvt_Rs[sort_idx]
    pvt_Bo = pvt_Bo[sort_idx]
    pvt_Bg = pvt_Bg[sort_idx]

    if metric:
        pvt_p = pvt_p * BAR_TO_PSI
        pvt_Rs = pvt_Rs * SM3_PER_SM3_TO_SCF_PER_STB
        # Bo: rm3/sm3 = rb/stb numerically (ratio of like units). No conversion.
        pvt_Bg = pvt_Bg / CUFTperBBL  # rm3/sm3 -> rb/scf
    return pvt_p, pvt_Rs, pvt_Bo, pvt_Bg


def _resolve_pvt(p_field, degf_field, pb_field, rsb_field, api, sg_g, sg_sp,
                 rsmethod, bomethod, zmethod, cmethod, pvt_arrays=None):
    """Oil/gas PVT at each survey pressure, in oilfield units (Bg as rb/scf).

    Interpolates the caller's table when pvt_arrays is supplied, otherwise walks
    the correlations one pressure at a time (the oil module is scalar-only).
    Returns (Rs, Bo, Bg, pb, rsb) with pb/rsb resolved from whichever source was
    used, since the caller needs them for the saturated/undersaturated split.
    """
    n = len(p_field)
    if pvt_arrays is not None:
        pvt_p, pvt_Rs, pvt_Bo, pvt_Bg = pvt_arrays
        if np.any(p_field < pvt_p[0] - 1e-6) or np.any(p_field > pvt_p[-1] + 1e-6):
            raise ValueError("Survey pressures must fall within pvt_table pressure range")
        Rs_arr = np.interp(p_field, pvt_p, pvt_Rs)
        Bo_arr = np.interp(p_field, pvt_p, pvt_Bo)
        Bg_arr = np.interp(p_field, pvt_p, pvt_Bg)
        # Infer pb/rsb from the table if not specified
        if pb_field <= 0:
            pb_field = pvt_p[np.argmax(pvt_Rs)]
        if rsb_field <= 0:
            rsb_field = np.max(pvt_Rs)
        return Rs_arr, Bo_arr, Bg_arr, pb_field, rsb_field

    sg_o = 141.5 / (api + 131.5)
    if pb_field <= 0 and rsb_field <= 0:
        raise ValueError("At least one of pb or rsb must be specified")
    if pb_field <= 0:
        pb_field = oil.oil_pbub(api, degf_field, rsb_field, sg_g=sg_g, sg_sp=sg_sp)
    if rsb_field <= 0:
        # Inverse of the oil_pbub (VALMC) relation above, so Pb and Rsb agree
        rsb_field = oil.oil_rs_bub(api, degf_field, pb_field, sg_g=sg_g, sg_sp=sg_sp,
                                    rsmethod='VALMC')

    Rs_arr = np.zeros(n)
    Bo_arr = np.zeros(n)
    Bg_arr = np.zeros(n)
    for i, pi in enumerate(p_field):
        rs_i = oil.oil_rs(api, degf_field, sg_sp, pi,
                          pb=pb_field, rsb=rsb_field, rsmethod=rsmethod)
        Rs_arr[i] = rs_i

        bo_i = oil.oil_bo(pi, pb_field, degf_field, rs_i, rsb_field, sg_o,
                          sg_g=sg_g, sg_sp=sg_sp, bomethod=bomethod)
        Bo_arr[i] = bo_i

        # Bg in rcf/scf, convert to rb/scf
        bg_i = gas.gas_bg(pi, sg_sp, degf_field, zmethod=zmethod, cmethod=cmethod) / CUFTperBBL
        Bg_arr[i] = bg_i

    return Rs_arr, Bo_arr, Bg_arr, pb_field, rsb_field


def oil_matbal(p, Np, degf, api=0, sg_sp=0, sg_g=0, pb=0, rsb=0,
               Rp=None, Wp=None, Wi=None, Gi=None, We=None,
               Bw=1.0, m=0, cf=0, sw_i=0, cw=0,
               rsmethod='VELAR', bomethod='MCAIN',
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None, regress=None):
    """Havlena-Odeh oil material balance for OOIP estimation.

    Parameters
    ----------
    p : array-like
        Reservoir pressures at each survey (psia | barsa). First value is
        initial reservoir pressure.
    Np : array-like
        Cumulative oil production (STB | sm3) at each pressure step. Same
        length as p. First entry is typically 0.
    degf : float
        Reservoir temperature (deg F | deg C).
    api : float
        Stock tank oil API gravity. Required when pvt_table is not provided.
    sg_sp : float
        Separator gas specific gravity. Required when pvt_table is not provided.
    sg_g : float
        Weighted average surface gas specific gravity (default 0 = use sg_sp).
    pb : float
        Bubble point pressure (psia | barsa). If 0, calculated from rsb.
    rsb : float
        Solution GOR at bubble point (scf/stb | sm3/sm3). If 0, calculated from pb.
    Rp : array-like, optional
        Cumulative producing GOR (scf/stb | sm3/sm3) at each step. Default = Rs.
    Wp : array-like, optional
        Cumulative water production (STB | sm3). Default all zeros.
    Wi : array-like, optional
        Cumulative water injection (STB | sm3). Default all zeros.
    Gi : array-like, optional
        Cumulative gas injection (scf | sm3). Default all zeros.
    We : array-like, optional
        Cumulative aquifer water influx in reservoir volume (rb | rm3). Same
        length as p. When provided, the influx is subtracted from underground
        withdrawal before estimating OOIP (Havlena-Odeh ``F - We = N·Et``), and a
        Water Drive Index ('WDI') is added to ``drive_indices``. Default all zeros.
    Bw : float
        Water FVF (rb/stb | rm3/sm3, default 1.0).
    m : float
        Gas cap size ratio (initial gas cap / initial oil zone PV, default 0).
    cf : float
        Formation compressibility (1/psi | 1/bar, default 0).
    sw_i : float
        Initial water saturation (fraction, default 0).
    cw : float
        Water compressibility (1/psi | 1/bar, default 0).
    rsmethod : str or rs_method
        Solution GOR method (default 'VELAR').
    bomethod : str or bo_method
        Oil FVF method (default 'MCAIN').
    zmethod : str or z_method
        Z-factor method for Bg (default 'DAK').
    cmethod : str or c_method
        Critical property method for Bg (default 'PMC').
    metric : bool
        If True, inputs/outputs in Eclipse METRIC units (default False).
    pvt_table : dict, optional
        Tabulated oil PVT: ``{'p': [...], 'Rs': [...], 'Bo': [...], 'Bg': [...]}``.
        Units follow the metric flag. When provided, api and sg_sp are not required.
    regress : dict, optional
        Parameters to regress with bounds: ``{'m': (0, 2), 'cf': (1e-6, 10e-6)}``.
        Allowed keys: 'm', 'cf', 'cw', 'sw_i'. Least squares on F - We - N Et with
        N the slope for each trial set. Warns when a parameter ends on a bound or
        is too correlated with N to be determined (see ``regress_report``).

    OOIP is the least-squares slope of F - We = N Et through the origin
    (Havlena and Odeh 1963).

    Returns
    -------
    OilMatbalResult
        With ``regressed`` and ``regress_report`` populated when regress is used.
    """
    p = np.asarray(p, dtype=float)
    Np = np.asarray(Np, dtype=float)
    n = len(p)

    if np.any(p <= 0):
        raise ValueError("Pressures must be positive")
    if len(Np) != n:
        raise ValueError(f"p and Np must have same length, got {len(p)} and {len(Np)}")
    if n < 2:
        raise ValueError("Need at least 2 pressure/production data points")

    _validate_regress(regress)

    # Validate PVT source
    if pvt_table is not None:
        pvt_p, pvt_Rs, pvt_Bo, pvt_Bg = _prepare_pvt_table(pvt_table, metric)
    elif api <= 0 or sg_sp <= 0:
        raise ValueError("api and sg_sp are required when pvt_table is not provided")

    # Convert metric to oilfield for internal calculations
    if metric:
        p_field = p * BAR_TO_PSI
        degf_field = degc_to_degf(degf)
        pb_field = pb * BAR_TO_PSI if pb > 0 else 0
        rsb_field = rsb * SM3_PER_SM3_TO_SCF_PER_STB if rsb > 0 else 0
        # Convert compressibilities from 1/bar to 1/psi
        cw = cw / BAR_TO_PSI   # 1/bar → 1/psi
        cf = cf / BAR_TO_PSI   # 1/bar → 1/psi
    else:
        p_field = p
        degf_field = degf
        pb_field = pb
        rsb_field = rsb

    # Default arrays
    if Rp is None:
        Rp_arr = None  # Will be set to Rs at each step
    else:
        Rp_arr = np.asarray(Rp, dtype=float)
        if metric:
            Rp_arr = Rp_arr * SM3_PER_SM3_TO_SCF_PER_STB

    Wp_arr = np.zeros(n) if Wp is None else np.asarray(Wp, dtype=float)
    Wi_arr = np.zeros(n) if Wi is None else np.asarray(Wi, dtype=float)
    Gi_arr = np.zeros(n) if Gi is None else np.asarray(Gi, dtype=float)
    if metric and Gi is not None:
        # Internal Bg is rb/scf (rm3/sm3 divided by CUFTperBBL), so the Gi*Bg
        # term needs Gi scaled by the same factor to recover a consistent
        # reservoir volume. SM3_PER_SM3_TO_SCF_PER_STB == CUFTperBBL numerically.
        Gi_arr = Gi_arr * SM3_PER_SM3_TO_SCF_PER_STB
    # Aquifer influx is already a reservoir volume (rb | rm3) matching F — no
    # unit conversion is needed in either metric or oilfield mode.
    We_arr = np.zeros(n) if We is None else np.asarray(We, dtype=float)
    for _name, _arr in (('Rp', Rp_arr), ('Wp', Wp_arr), ('Wi', Wi_arr),
                        ('Gi', Gi_arr), ('We', We_arr)):
        if _arr is not None and len(_arr) != n:
            raise ValueError(f"{_name} must have same length as p, got {len(_arr)} and {n}")

    pvt_arrays = None if pvt_table is None else (pvt_p, pvt_Rs, pvt_Bo, pvt_Bg)
    Rs_arr, Bo_arr, Bg_arr, pb_field, rsb_field = _resolve_pvt(
        p_field, degf_field, pb_field, rsb_field, api, sg_g, sg_sp,
        rsmethod, bomethod, zmethod, cmethod, pvt_arrays)

    # Initial conditions
    Boi = Bo_arr[0]
    Rsi = Rs_arr[0]
    Bgi = Bg_arr[0]

    # If Rp not provided, use Rs (no free gas produced)
    if Rp_arr is None:
        Rp_arr = Rs_arr.copy()

    # Precompute F, Eo, Eg (these do NOT depend on regressable params)
    F = (Np * (Bo_arr + (Rp_arr - Rs_arr) * Bg_arr)
         + (Wp_arr - Wi_arr) * Bw
         - Gi_arr * Bg_arr)

    # Net withdrawal driving oil/gas/rock expansion: aquifer influx supplies part
    # of the withdrawn volume, so it is removed before estimating OOIP
    # (Havlena-Odeh: F - We = N * [Eo + m*Eg + (1+m)*Efw]).
    F_net = F - We_arr

    Eo = np.where(p_field >= pb_field,
                  Bo_arr - Boi,
                  (Bo_arr - Boi) + (Rsi - Rs_arr) * Bg_arr)

    Eg = np.where(Bgi > 0, Boi * (Bg_arr / max(Bgi, 1e-30) - 1.0), 0.0)
    dp = p_field[0] - p_field

    regressed_result = report = None
    if regress is not None:
        fitted, report = _regress_oil(regress, {'m': m, 'cf': cf, 'cw': cw, 'sw_i': sw_i},
                                      F_net, Eo, Eg, dp, Boi, metric)
        m, cf, cw, sw_i = (fitted.get(k, v) for k, v in
                           (('m', m), ('cf', cf), ('cw', cw), ('sw_i', sw_i)))
        # Regressed compressibilities back in the caller's units (1/psi -> 1/bar)
        regressed_result = {k: float(v * BAR_TO_PSI if metric and k in ('cw', 'cf') else v)
                            for k, v in fitted.items()}

    Et, Efw = _havlena_odeh_et(Eo, Eg, dp, Boi, m, cf, cw, sw_i)
    ooip = _ooip_slope(F_net, Et)
    if ooip is None:
        raise RuntimeError("Cannot compute OOIP — all expansion terms are zero")

    # Drive indices as fractions of the total withdrawal F, from the fitted N: their
    # sum is 1 only where N Et + We reproduces F (Pletcher 2002).
    has_f = np.abs(F) > 1e-30
    F_safe = np.where(has_f, F, 1.0)
    ddi, sdi, cdi, wdi = (np.where(has_f, term / F_safe, 0.0) for term in
                          (ooip * Eo, ooip * m * Eg, ooip * (1.0 + m) * Efw, We_arr))
    has_et = np.abs(Et) > 1e-30
    ooip_points = np.where(has_et, F_net / np.where(has_et, Et, 1.0), np.nan)

    return OilMatbalResult(
        ooip=ooip,
        F=F,
        Eo=Eo,
        Eg=Eg,
        Efw=Efw,
        drive_indices={'DDI': ddi, 'SDI': sdi, 'CDI': cdi, 'WDI': wdi},
        p=p,
        pvt={'Rs': Rs_arr, 'Bo': Bo_arr, 'Bg': Bg_arr},
        regressed=regressed_result,
        ooip_points=ooip_points,
        regress_report=report,
    )
