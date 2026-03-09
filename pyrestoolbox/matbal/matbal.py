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

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List

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
        Drive index fractions at each step: 'DDI' (depletion), 'SDI' (segregation/gas cap),
        'CDI' (compaction/water). Each maps to np.ndarray.
    p : np.ndarray
        Pressure array used.
    pvt : dict
        PVT properties at each step: 'Rs', 'Bo', 'Bg' arrays.
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


def gas_matbal(p, Gp, degf, sg=0.65, co2=0, h2s=0, n2=0, h2=0,
               Wp=None, Bw=1.0, We=None,
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None):
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
    Wp : array-like, optional
        Cumulative water production (STB | sm3). Same length as p.
    Bw : float
        Water FVF (rb/stb | rm3/sm3, default 1.0). Used when Wp provided.
    We : array-like, optional
        Cumulative water influx (rcf | rm3). Same length as p. When
        provided, Havlena-Odeh regression is used for OGIP instead of P/Z.
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
            gas.gas_z(pi, sg, degf, co2=co2, h2s=h2s, n2=n2, h2=h2,
                      zmethod=zmethod, cmethod=cmethod, metric=metric)
            for pi in p
        ])

        # Bg at each pressure (rcf/scf — gas_bg handles metric internally)
        bg = np.array([
            gas.gas_bg(pi, sg, degf, co2=co2, h2s=h2s, n2=n2, h2=h2,
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
    )


def oil_matbal(p, Np, degf, api=0, sg_sp=0, sg_g=0, pb=0, rsb=0,
               Rp=None, Wp=None, Wi=None, Gi=None,
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
        Allowed keys: 'm', 'cf', 'cw', 'sw_i'. Optimizes to minimize the
        coefficient of variation of OOIP estimates across time steps.

    Returns
    -------
    OilMatbalResult
        With ``regressed`` dict populated when regress is used.
    """
    p = np.asarray(p, dtype=float)
    Np = np.asarray(Np, dtype=float)
    n = len(p)

    if len(Np) != n:
        raise ValueError(f"p and Np must have same length, got {len(p)} and {len(Np)}")
    if n < 2:
        raise ValueError("Need at least 2 pressure/production data points")

    # Validate regress parameter
    _ALLOWED_REGRESS = {'m', 'cf', 'cw', 'sw_i'}
    if regress is not None:
        for key, val in regress.items():
            if key not in _ALLOWED_REGRESS:
                raise ValueError(f"Invalid regress key '{key}'. Allowed: {sorted(_ALLOWED_REGRESS)}")
            if not isinstance(val, (tuple, list)) or len(val) != 2:
                raise ValueError(f"regress['{key}'] must be a (lower, upper) tuple, got {val}")
            if val[0] >= val[1]:
                raise ValueError(f"regress['{key}'] bounds must have lower < upper, got {val}")

    # Validate PVT source
    if pvt_table is not None:
        # Validate pvt_table keys
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
        # Sort by pressure
        sort_idx = np.argsort(pvt_p)
        pvt_p = pvt_p[sort_idx]
        pvt_Rs = pvt_Rs[sort_idx]
        pvt_Bo = pvt_Bo[sort_idx]
        pvt_Bg = pvt_Bg[sort_idx]
        # Metric conversion
        if metric:
            pvt_p = pvt_p * BAR_TO_PSI
            pvt_Rs = pvt_Rs * SM3_PER_SM3_TO_SCF_PER_STB
            # Bo: rm3/sm3 = rb/stb numerically (ratio of like units). No conversion.
            pvt_Bg = pvt_Bg / CUFTperBBL  # rm3/sm3 → rb/scf
    elif api <= 0 or sg_sp <= 0:
        raise ValueError("api and sg_sp are required when pvt_table is not provided")

    # Convert metric to oilfield for internal calculations
    if metric:
        p_field = p * BAR_TO_PSI
        degf_field = degc_to_degf(degf)
        pb_field = pb * BAR_TO_PSI if pb > 0 else 0
        rsb_field = rsb * SM3_PER_SM3_TO_SCF_PER_STB if rsb > 0 else 0
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

    if pvt_table is not None:
        # Check survey pressures within table range
        if np.any(p_field < pvt_p[0] - 1e-6) or np.any(p_field > pvt_p[-1] + 1e-6):
            raise ValueError("Survey pressures must fall within pvt_table pressure range")
        # Interpolate PVT from table
        Rs_arr = np.interp(p_field, pvt_p, pvt_Rs)
        Bo_arr = np.interp(p_field, pvt_p, pvt_Bo)
        Bg_arr = np.interp(p_field, pvt_p, pvt_Bg)
        # Infer pb/rsb from table if not specified
        if pb_field <= 0:
            pb_field = pvt_p[np.argmax(pvt_Rs)]
        if rsb_field <= 0:
            rsb_field = np.max(pvt_Rs)
    else:
        sg_o = 141.5 / (api + 131.5)
        # Resolve pb and rsb
        if pb_field <= 0 and rsb_field <= 0:
            raise ValueError("At least one of pb or rsb must be specified")
        if pb_field <= 0:
            pb_field = oil.oil_pbub(api, degf_field, rsb_field, sg_g=sg_g, sg_sp=sg_sp)
        if rsb_field <= 0:
            rsb_field = oil.oil_rs_bub(api, degf_field, pb_field, sg_g=sg_g, sg_sp=sg_sp,
                                        rsmethod=rsmethod)
        # Compute PVT at each pressure (oil module is scalar-only)
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

    # Initial conditions
    Boi = Bo_arr[0]
    Rsi = Rs_arr[0]
    Bgi = Bg_arr[0]

    # If Rp not provided, use Rs (no free gas produced)
    if Rp_arr is None:
        Rp_arr = Rs_arr.copy()

    # Precompute F, Eo, Eg (these do NOT depend on regressable params)
    F = np.zeros(n)
    Eo = np.zeros(n)
    Eg = np.zeros(n)

    for i in range(n):
        # Underground withdrawal
        F[i] = (Np[i] * (Bo_arr[i] + (Rp_arr[i] - Rs_arr[i]) * Bg_arr[i])
                + (Wp_arr[i] - Wi_arr[i]) * Bw
                - Gi_arr[i] * Bg_arr[i])

        # Oil expansion + dissolved gas
        if p_field[i] >= pb_field:
            Eo[i] = Bo_arr[i] - Boi
        else:
            Eo[i] = (Bo_arr[i] - Boi) + (Rsi - Rs_arr[i]) * Bg_arr[i]

        # Gas cap expansion
        if Bgi > 0:
            Eg[i] = Boi * (Bg_arr[i] / Bgi - 1.0)
        else:
            Eg[i] = 0.0

    # Helper to compute Efw and OOIP for given (m_t, cf_t, cw_t, sw_t)
    def _compute_result(m_t, cf_t, cw_t, sw_t):
        Efw = np.zeros(n)
        for i in range(n):
            dp = p_field[0] - p_field[i]
            if (1.0 - sw_t) > 0:
                Efw[i] = Boi * (cw_t * sw_t + cf_t) / (1.0 - sw_t) * dp
            else:
                Efw[i] = 0.0

        denom = Eo + m_t * Eg + (1.0 + m_t) * Efw
        valid = np.abs(denom) > 1e-30
        if not np.any(valid):
            return None, Efw, denom, valid
        N_estimates = F[valid] / denom[valid]
        ooip = float(np.mean(N_estimates))
        return ooip, Efw, denom, valid

    # Regression
    regressed_result = None
    if regress is not None:
        from scipy.optimize import minimize

        param_names = list(regress.keys())
        bounds = [regress[k] for k in param_names]
        base_vals = {'m': m, 'cf': cf, 'cw': cw, 'sw_i': sw_i}

        # Initial guess: use passed value if within bounds, else midpoint
        x0 = []
        for k, (lo, hi) in zip(param_names, bounds):
            v = base_vals[k]
            if lo <= v <= hi:
                x0.append(v)
            else:
                x0.append((lo + hi) / 2.0)

        def objective(x):
            vals = dict(base_vals)
            for k, v in zip(param_names, x):
                vals[k] = v
            Efw_t = np.zeros(n)
            for i in range(n):
                dp = p_field[0] - p_field[i]
                if (1.0 - vals['sw_i']) > 0:
                    Efw_t[i] = Boi * (vals['cw'] * vals['sw_i'] + vals['cf']) / (1.0 - vals['sw_i']) * dp
            denom_t = Eo + vals['m'] * Eg + (1.0 + vals['m']) * Efw_t
            valid_t = np.abs(denom_t) > 1e-30
            if not np.any(valid_t):
                return 1e10
            N_est = F[valid_t] / denom_t[valid_t]
            mean_N = np.mean(N_est)
            if abs(mean_N) < 1e-30:
                return 1e10
            return float(np.std(N_est) / abs(mean_N))

        opt_result = minimize(objective, x0, method='L-BFGS-B', bounds=bounds)
        # Apply optimized values
        for k, v in zip(param_names, opt_result.x):
            base_vals[k] = v
        m = base_vals['m']
        cf = base_vals['cf']
        cw = base_vals['cw']
        sw_i = base_vals['sw_i']
        regressed_result = {k: float(v) for k, v in zip(param_names, opt_result.x)}

    # Final computation with (possibly regressed) params
    ooip, Efw, denom, valid = _compute_result(m, cf, cw, sw_i)
    if ooip is None:
        raise RuntimeError("Cannot compute OOIP — all expansion terms are zero")

    # Drive indices at each valid point
    ddi = np.zeros(n)
    sdi = np.zeros(n)
    cdi = np.zeros(n)
    for i in range(n):
        if abs(denom[i]) > 1e-30 and abs(F[i]) > 1e-30:
            N_i = F[i] / denom[i]
            ddi[i] = N_i * Eo[i] / F[i]
            sdi[i] = N_i * m * Eg[i] / F[i]
            cdi[i] = N_i * (1.0 + m) * Efw[i] / F[i]

    return OilMatbalResult(
        ooip=ooip,
        F=F,
        Eo=Eo,
        Eg=Eg,
        Efw=Efw,
        drive_indices={'DDI': ddi, 'SDI': sdi, 'CDI': cdi},
        p=p,
        pvt={'Rs': Rs_arr, 'Bo': Bo_arr, 'Bg': Bg_arr},
        regressed=regressed_result,
    )
