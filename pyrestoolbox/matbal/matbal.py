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
from pyrestoolbox.constants import CUFTperBBL
from pyrestoolbox.constants import BAR_TO_PSI, degc_to_degf, PSI_TO_BAR
from pyrestoolbox.constants import SCF_PER_STB_TO_SM3_PER_SM3, SM3_PER_SM3_TO_SCF_PER_STB
from pyrestoolbox.classes import z_method, c_method, rs_method, bo_method


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
    """
    ogip: float
    pz: np.ndarray
    gp: np.ndarray
    slope: float
    intercept: float
    r_squared: float
    p_initial: float
    z_initial: float


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


def gas_matbal(p, Gp, degf, sg=0.65, co2=0, h2s=0, n2=0, h2=0,
               zmethod='DAK', cmethod='PMC', metric=False):
    """P/Z gas material balance for OGIP estimation.

    Performs linear regression of P/Z vs cumulative gas production to
    determine original gas in place (OGIP = -intercept/slope).

    Parameters
    ----------
    p : array-like
        Reservoir pressures at each survey (psia | barsa). First value is
        initial pressure (Gp=0 at that point, or the first Gp entry).
    Gp : array-like
        Cumulative gas production at each pressure survey. Same length as p.
        Units are user-defined (e.g. Bscf, MMscf) — OGIP will be in the same units.
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
    zmethod : str or z_method
        Z-factor method (default 'DAK').
    cmethod : str or c_method
        Critical property method (default 'PMC').
    metric : bool
        If True, p in barsa and degf in deg C (default False).

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

    # Convert metric inputs for gas_z (which handles metric internally)
    # gas_z signature: gas_z(p, sg, degf, ...)
    z = np.array([
        gas.gas_z(pi, sg, degf, co2=co2, h2s=h2s, n2=n2, h2=h2,
                  zmethod=zmethod, cmethod=cmethod, metric=metric)
        for pi in p
    ])

    pz = p / z

    # Linear regression: P/Z = intercept + slope * Gp
    coeffs = np.polyfit(Gp, pz, 1)
    slope = coeffs[0]
    intercept = coeffs[1]

    # OGIP = -intercept / slope (where P/Z = 0)
    if abs(slope) < 1e-30:
        raise RuntimeError("P/Z vs Gp regression has near-zero slope — cannot determine OGIP")
    ogip = -intercept / slope

    # R-squared
    pz_pred = np.polyval(coeffs, Gp)
    ss_res = np.sum((pz - pz_pred) ** 2)
    ss_tot = np.sum((pz - np.mean(pz)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return GasMatbalResult(
        ogip=float(ogip),
        pz=pz,
        gp=Gp,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
        p_initial=float(p[0]),
        z_initial=float(z[0]),
    )


def oil_matbal(p, Np, degf, api, sg_sp, sg_g=0, pb=0, rsb=0,
               Rp=None, Wp=None, Wi=None, Gi=None,
               Bw=1.0, m=0, cf=0, sw_i=0, cw=0,
               rsmethod='VELAR', bomethod='MCAIN',
               zmethod='DAK', cmethod='PMC', metric=False):
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
        Stock tank oil API gravity.
    sg_sp : float
        Separator gas specific gravity.
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

    Returns
    -------
    OilMatbalResult
    """
    p = np.asarray(p, dtype=float)
    Np = np.asarray(Np, dtype=float)
    n = len(p)

    if len(Np) != n:
        raise ValueError(f"p and Np must have same length, got {len(p)} and {len(Np)}")
    if n < 2:
        raise ValueError("Need at least 2 pressure/production data points")

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

    sg_o = 141.5 / (api + 131.5)

    # Resolve pb and rsb
    if pb_field <= 0 and rsb_field <= 0:
        raise ValueError("At least one of pb or rsb must be specified")
    if pb_field <= 0:
        pb_field = oil.oil_pbub(api, degf_field, rsb_field, sg_g=sg_g, sg_sp=sg_sp)
    if rsb_field <= 0:
        rsb_field = oil.oil_rs_bub(api, degf_field, pb_field, sg_g=sg_g, sg_sp=sg_sp,
                                    rsmethod=rsmethod)

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
        # gas_bg signature: gas_bg(p, sg, degf, ...)
        bg_i = gas.gas_bg(pi, sg_sp, degf_field, zmethod=zmethod, cmethod=cmethod) / CUFTperBBL
        Bg_arr[i] = bg_i

    # Initial conditions
    Boi = Bo_arr[0]
    Rsi = Rs_arr[0]
    Bgi = Bg_arr[0]

    # If Rp not provided, use Rs (no free gas produced)
    if Rp_arr is None:
        Rp_arr = Rs_arr.copy()

    # Havlena-Odeh material balance
    F = np.zeros(n)
    Eo = np.zeros(n)
    Eg = np.zeros(n)
    Efw = np.zeros(n)

    for i in range(n):
        dp = p_field[0] - p_field[i]

        # Underground withdrawal
        F[i] = (Np[i] * (Bo_arr[i] + (Rp_arr[i] - Rs_arr[i]) * Bg_arr[i])
                + (Wp_arr[i] - Wi_arr[i]) * Bw
                - Gi_arr[i] * Bg_arr[i])

        # Oil expansion + dissolved gas
        if p_field[i] >= pb_field:
            # Above bubble point: no free gas, Rs = Rsi
            Eo[i] = Bo_arr[i] - Boi
        else:
            Eo[i] = (Bo_arr[i] - Boi) + (Rsi - Rs_arr[i]) * Bg_arr[i]

        # Gas cap expansion
        if Bgi > 0:
            Eg[i] = Boi * (Bg_arr[i] / Bgi - 1.0)
        else:
            Eg[i] = 0.0

        # Formation and water compressibility
        if (1.0 - sw_i) > 0:
            Efw[i] = Boi * (cw * sw_i + cf) / (1.0 - sw_i) * dp
        else:
            Efw[i] = 0.0

    # Compute OOIP: N = F / (Eo + m*Eg + (1+m)*Efw)
    denom = Eo + m * Eg + (1.0 + m) * Efw

    # Use points where denom > 0 (skip initial point where everything is zero)
    valid = np.abs(denom) > 1e-30
    if not np.any(valid):
        raise RuntimeError("Cannot compute OOIP — all expansion terms are zero")

    # Weighted average OOIP (F/denom at each valid point)
    N_estimates = F[valid] / denom[valid]
    ooip = float(np.mean(N_estimates))

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
    )
