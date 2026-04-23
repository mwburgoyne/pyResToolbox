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

Gas PVT and flow rate calculations.

Functions
---------
gas_z           Z-factor (DAK, HY, WYW, or BNS methods)
gas_ug          Gas viscosity (cP)
gas_bg          Gas formation volume factor (rcf/scf)
gas_cg          Gas compressibility (1/psi)
gas_den         Gas density (lb/cuft)
gas_sg          Gas specific gravity from composition
gas_tc_pc       Pseudocritical temperature & pressure
gas_dmp         Gas pseudopressure difference
gas_ponz2p      Convert P/Z to pressure
gas_grad2sg     Invert gas gradient to SG
gas_fws_sg      Fresh-water-saturated gas SG
gas_water_content  Equilibrium water content of gas
gas_rate_radial Radial gas flow rate (Mscf/d)
gas_rate_linear Linear gas flow rate (Mscf/d)
darcy_gas       Darcy gas rate from pseudopressure difference
gas_hydrate     Gas hydrate formation prediction and inhibitor calculations
gas_hvf_beta    Forchheimer β via Firoozabadi-Katz / Jones / Tek-Coats-Katz
gas_non_darcy_skin          Rate-dependent (HVF) skin S_hvf = D·q_g
gas_partial_penetration_skin  Streltsova-Adams partial-penetration pseudoskin

Classes
-------
GasPVT          Convenience wrapper storing gas composition & method choices
HydrateResult   Dataclass returned by gas_hydrate()
"""

__all__ = [
    'gas_z', 'gas_ug', 'gas_bg', 'gas_cg', 'gas_den', 'gas_sg',
    'gas_tc_pc', 'gas_dmp', 'gas_ponz2p', 'gas_grad2sg', 'gas_fws_sg',
    'gas_water_content', 'gas_rate_radial', 'gas_rate_linear', 'darcy_gas',
    'gas_hydrate',
    'gas_hvf_beta', 'gas_non_darcy_skin', 'gas_partial_penetration_skin',
    'GasPVT', 'HydrateResult',
    # Enum classes re-exported for convenience (gas.z_method.DAK, etc.)
    'z_method', 'c_method', 'hyd_method', 'inhibitor',
]

import math
import warnings
import numpy as np
import numpy.typing as npt
from typing import Tuple
from dataclasses import dataclass

from pyrestoolbox.classes import z_method, c_method, hyd_method, inhibitor
from pyrestoolbox.shared_fns import convert_to_numpy, process_output, check_2_inputs, bisect_solve, validate_pe_inputs
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.constants import (R, psc, tsc, degF2R, scf_per_mol, CUFTperBBL, WDEN, MW_CO2, MW_H2S, MW_N2, MW_AIR, MW_H2,
    BAR_TO_PSI, PSI_TO_BAR, degc_to_degf, degf_to_degc,
    M_TO_FT, SQM_TO_SQFT,
    LBCUFT_TO_KGM3, INVPSI_TO_INVBAR,
    PSI2CP_TO_BAR2CP, BARM_TO_PSIFT,
    MSCF_TO_SM3, STB_PER_MMSCF_TO_SM3_PER_SM3,
    SM3_PER_SM3_TO_STB_PER_MMSCF,
    D_PER_SM3_TO_D_PER_MSCF,
    LB_PER_MMSCF_TO_KG_PER_SM3, GAL_PER_MMSCF_TO_L_PER_SM3)

# Precomputed Gauss-Legendre nodes/weights for pseudopressure integration
_GL7_NODES, _GL7_WEIGHTS = np.polynomial.legendre.leggauss(7)
_GL10_NODES, _GL10_WEIGHTS = np.polynomial.legendre.leggauss(10)

# gas_grad2sg bisection bounds. Lower = pure H2 SG (physical minimum for
# H2-blend support); upper = 3.0 covers pure CO2 (SG ~1.53) with margin.
_GRAD2SG_SG_LO = MW_H2 / MW_AIR
_GRAD2SG_SG_HI = 3.0

# =============================================================================
# Correlation coefficients — named constants with paper citations
# =============================================================================

# --- Piper, McCain & Corredor (1999) ---
# Eqs 2.4-2.6, 'Petroleum Reservoir Fluid Property Correlations', McCain et al.
_PMC_ALPHA = np.array([0.11582, -0.4582, -0.90348, -0.66026, 0.70729, -0.099397])
_PMC_BETA = np.array([3.8216, -0.06534, -0.42113, -0.91249, 17.438, -3.2191])
_PMC_TCI = np.array([0, 672.35, 547.58, 239.26])   # H2S, CO2, N2 critical T (deg R)
_PMC_PCI = np.array([0, 1306.0, 1071.0, 507.5])    # H2S, CO2, N2 critical P (psia)

# --- Sutton (1985) with Wichert & Aziz (1972) corrections ---
# Eqs 3.47a-b, 3.52-3.54, McCain et al.
_SUT_TPC = (169.2, 349.5, -74.0)           # Tpc = a + b*sg + c*sg^2 (deg R)
_SUT_PPC = (756.8, -131.0, -3.6)           # Ppc = a + b*sg + c*sg^2 (psia)
_SUT_WA = (120.0, 0.9, 1.6, 15.0, 0.5, 4.0)  # Wichert-Aziz epsilon coefficients
# Component critical properties for SUT/PMC mixing (N2, CO2, H2S) — deg R and psia
_IMPUR_TC = (239.26, 547.58, 672.35)       # N2, CO2, H2S critical T (deg R)
_IMPUR_PC = (507.5, 1071.0, 1306.0)        # N2, CO2, H2S critical P (psia)

# --- BNS pseudo-critical fits (Burgoyne, 2025) ---
_BNS_TC_AG = (2695.14765, 274.341701, 343.008)   # Associated gas Tc Michaelis-Menten
_BNS_TC_GC = (1098.10948, 101.529237, 343.008)   # Gas condensate Tc Michaelis-Menten
_BNS_VC_INTERCEPT = 5.518525872412144             # Vc/Zc intercept for pc_fn
_BNS_VC_SLOPE_AG = 0.177497835                    # Associated gas Vc slope
_BNS_VC_SLOPE_GC = 0.170931432                    # Gas condensate Vc slope
_BNS_CH4_MW = 16.0425                             # Methane MW for BNS fits
_BNS_SG_METHANE = 0.553779772                     # Methane SG lower bound
_BNS_VCVIS_HC = (0.0576710, 1.44383)              # HC Vc linear fit: slope, intercept

# --- Dranchuk & Abou-Kassem (1975) ---
# Eqs 2.7-2.8, 'Petroleum Reservoir Fluid Property Correlations', McCain et al.
_DAK_A = (0.3265, -1.0700, -0.5339, 0.01569, -0.05165,
          0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210)
_DAK_LEADING = 0.27  # Leading coefficient in DAK reduced density equation

# --- Hall & Yarborough (1973) ---
_HY_COEFFS = (0.06125, -1.2, 14.76, -9.76, 4.58, 90.7, -242.2, 42.4, 2.18, 2.82)

# --- Wang, Ye & Wu (2021) ---
# doi:10.1016/j.egyr.2021.11.029
_WYW_A = np.array([0, 256.41675, 7.18202, -178.5725, 182.98704, -40.74427,
                    2.24427, 47.44825, 5.2852, -0.14914, 271.50446, 16.2694,
                    -121.51728, 167.71477, -81.73093, 20.36191, -2.1177,
                    124.64444, -6.74331, 0.20897, -0.00314])

# --- Lee, Gonzalez & Eakin (1966) ---
# Eqs 2.14-2.17, 'Petroleum Reservoir Fluid Property Correlations', McCain et al.
_LGE = (3.448, 986.4, 0.01009,    # b coefficients (Eq 2.16)
        2.447, 0.2224,              # c coefficients (Eq 2.17)
        9.379, 0.01607, 209.2, 19.26,  # a coefficients (Eq 2.15)
        0.0001)                     # output scale factor (Eq 2.14)

# --- Stiel-Thodos dilute gas viscosity ---
_ST_LOW = (34e-5, 0.94)              # Low-Tr: coeff, exponent
_ST_HIGH = (17.78e-5, 4.58, -1.67)   # High-Tr: coeff, slope, intercept
_ST_TR_THRESH = 1.5                   # Tr threshold between low/high

# --- Lohrenz-Bray-Clark (LBC) residual viscosity ---
_LBC_POLY = np.array([0.1023, 0.023364, 0.058533, -0.0392852, 0.00926279])
_LBC_OFFSET = 1e-4

# --- Darcy gas flow equation ---
_DARCY_GAS_CONST = 1422  # Field-unit conversion factor for gas pseudopressure flow

# --- Danesh water content correlation ---
# 'PVT and Phase Behaviour Of Petroleum Reservoir Fluids', Danesh
_DANESH_WC = (47484.0, 69.103501, -13064.76, -7.3037, 0.0000012856,  # vapour pressure term
              -3083.87, 6.69449,                                        # empirical second term
              0.00492, 0.00017672,                                      # salinity correction
              8.32, 42.0)                                               # unit conversion: lb/gal, gal/bbl

# --- Standing condensate MW correlation ---
_STANDING_COND_MW = (240.0, 2.22)  # MW = a - b * API

# --- Motiee (1991) hydrate formation temperature ---
# Hydrocarbon Processing 70, pp 98-99
_MOTIEE = (-283.24469, 78.99667, -5.352544, 349.473877, -150.854675, -27.604065)

# --- Towler & Mokhatab (2005) hydrate formation temperature ---
# Hydrocarbon Processing 84, pp 61-62
_TOWLER = (13.47, 34.27, -1.675, -20.35)

# --- Hydrate formation pressure search bounds ---
_HFP_P_LO = 14.696   # ~1 atm (psia)
_HFP_P_HI = 15000.0   # Upper search bound (psia)

# --- Non-Darcy / high-velocity-flow skin correlations ---
# "FK": simple log-log fit of the consolidated-rock β(k) data published by
# Firoozabadi & Katz (1979) "An Analysis of High-Velocity Gas Flow Through
# Porous Media," JPT Feb 1979, pp.211-216 (SPE-6827):
#   β[1/ft] = _FK_A * k^(-_FK_N), k in md
_FK_A = 2.172e10
_FK_N = 1.201

# Jones, S.C. (1987) "Using the Inertial Coefficient b to Characterize
# Heterogeneity in Reservoir Rock," SPE-16949:
#   β[1/ft] = _JONES_A * k^(-_JONES_N), k in md
_JONES_A = 6.15e10
_JONES_N = 1.55

# Tek, M.R., Coats, K.H., Katz, D.L. (1962) "The Effect of Turbulence on
# Flow of Natural Gas through Porous Reservoirs," JPT July 1962, pp.799-806:
#   β[1/ft] = _TCK_A / (k^_TCK_NK * φ^_TCK_NPHI), k in md, φ fraction
_TCK_A = 1.88e10
_TCK_NK = 1.47
_TCK_NPHI = 0.53

# Non-Darcy coefficient in field units -- Jones (1987) SPE-16949; also
# derivable from Odeh, Moreland & Schueler (1975) "Characterization of a
# Gas Well From One Flow-Test Sequence," JPT Dec 1975, pp.1501-1504.
#   D[day/MSCF] = _D_COATS_KATZ * β[1/ft] * γg * k[md] / (μg[cp] * h[ft] * rw[ft])
_D_COATS_KATZ = 2.222e-15

# Streltsova-Adams partial-penetration series: sum up to _SP_MAX_N terms.
# K_0 decays exponentially for large n·π·rD, but rD ≈ rw/h is typically ~1e-3
# for gas wells, so convergence can require tens of thousands of terms.
# Vectorised evaluation makes this sub-millisecond even at the cap.
_SP_MAX_N = 20000
_SP_CONV_TOL = 1.0e-6  # relative tail-block tolerance used to emit a warning

# Valid β-correlation tags
_BETA_METHODS = ('FK', 'JONES', 'TCK')

def _h2_method_override(h2, zmethod, cmethod):
    """Force BNS method when hydrogen is present."""
    if h2 > 0:
        return 'BNS', 'BNS'
    return zmethod, cmethod

def _is_bns_method(method):
    """True if method is BNS/BUR (Enum, string, or legacy alias)."""
    if hasattr(method, 'name'):
        return method.name in ('BNS', 'BUR')
    if isinstance(method, str):
        return method.upper() in ('BNS', 'BUR')
    return False

def _method_label(method):
    """Printable label for a method arg (Enum.name or string)."""
    if hasattr(method, 'name'):
        return method.name
    return str(method)

def _resolve_methods(zmethod, cmethod, h2=0):
    """Resolve z/c methods with BNS coupling and return validated Enums.

    Policy:
      1. h2 > 0 forces both methods to 'BNS' (documented auto-selection, no warning).
      2. If either method is BNS, force both to BNS and emit UserWarning naming
         the overruled counterpart. Non-BNS methods are not coupled.

    User-supplied tc/pc are always respected within the resulting method:
    for BNS they replace the hydrocarbon pseudo-component Tc/Pc (inert Tc/Pc
    remain BNS internal constants); for SUT/PMC they replace the mixture Tc/Pc.
    """
    if h2 > 0:
        zmethod, cmethod = 'BNS', 'BNS'
    else:
        z_is_bns = _is_bns_method(zmethod)
        c_is_bns = _is_bns_method(cmethod)
        if z_is_bns != c_is_bns:
            if z_is_bns:
                overruled = f"cmethod={_method_label(cmethod)!r}"
                cmethod = 'BNS'
            else:
                overruled = f"zmethod={_method_label(zmethod)!r}"
                zmethod = 'BNS'
            warnings.warn(
                f"BNS coupling: {overruled} overruled to 'BNS' because the "
                f"counterpart method is BNS. When any of zmethod/cmethod is BNS, "
                f"both are forced to BNS for thermodynamic consistency.",
                UserWarning,
                stacklevel=3,
            )
    return validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

def _metric_to_field_pvt(p, degf, tc, pc, metric):
    """Convert gas PVT inputs from metric (barsa, degC, K, barsa) to field (psia, degF, degR, psia)."""
    if not metric:
        return p, degf, tc, pc
    p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
    degf = degc_to_degf(degf)
    if tc > 0:
        tc = tc * 1.8  # K to deg R
    if pc > 0:
        pc = pc * BAR_TO_PSI
    return p, degf, tc, pc

# Optional Rust acceleration
from pyrestoolbox._accelerator import RUST_AVAILABLE, _rust_module

def _compute_delta_mp(pr, pwf, degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s):
    """Compute pseudopressure difference and flow direction for gas rate calculations.

    Handles scalar and array inputs for pr and pwf.
    Returns (direction, delta_mp) where direction encodes flow sign.
    """
    direction = 1
    if pr.size + pwf.size == 2:  # Single set of pressures
        if pr < pwf:
            direction = -1
        delta_mp = abs(
            gas_dmp(
                p1=pwf, p2=pr, degf=degf, sg=sg,
                zmethod=zmethod, cmethod=cmethod, tc=tc, pc=pc,
                n2=n2, co2=co2, h2s=h2s,
            )
        )
    else:
        if pr.size > 1:  # Multiple Pr's
            direction = 2 * np.array([p > pwf for p in pr]) - 1
            delta_mp = np.absolute(np.array([
                gas_dmp(p1=p, p2=pwf, degf=degf, sg=sg,
                        zmethod=zmethod, cmethod=cmethod, tc=tc, pc=pc,
                        n2=n2, co2=co2, h2s=h2s)
                for p in pr
            ]))
        else:  # Multiple BHFP's
            direction = 2 * np.array([pr > bhfp for bhfp in pwf]) - 1
            delta_mp = np.absolute(np.array([
                gas_dmp(p1=pr, p2=bhfp, degf=degf, sg=sg,
                        zmethod=zmethod, cmethod=cmethod, tc=tc, pc=pc,
                        n2=n2, co2=co2, h2s=h2s)
                for bhfp in pwf
            ]))
    return direction, delta_mp

def _prepare_gas_rate_inputs(degf, sg, co2, h2s, n2, h2, zmethod, cmethod, tc, pc):
    """Validate inputs, apply H2 auto-selection, resolve methods and critical properties."""
    validate_pe_inputs(degf=degf, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)
    tc, pc = gas_tc_pc(sg, co2, h2s, n2, h2, cmethod.name, tc, pc)
    return zmethod, cmethod, tc, pc

def gas_rate_radial(
    k: npt.ArrayLike,
    h: npt.ArrayLike,
    pr: npt.ArrayLike,
    pwf: npt.ArrayLike,
    r_w: float,
    r_ext: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    S: float = 0,
    D: float = 0,
    sg: float = 0.75,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    gas_pvt = None,
    metric: bool = False,
) -> np.ndarray:
    """ Returns gas rate for radial flow (mscf/day) using Darcy pseudo steady state equation & gas pseudopressure
        k: Permeability (mD)
        h: Net flow height (ft)
        pr: Reservoir pressure (psia)
        pwf: BHFP (psia)
        r_w: Wellbore Radius (ft)
        r_ext: External Reservoir Radius (ft)
        degf: Reservoir Temperature (deg F)
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                 defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculating critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
               Defaults to 'PMC'
                S: Skin. Defaults to zero if undefined
        D: Non Darcy Skin Factor (day/mscf). Defaults to zero if undefined
        sg: Gas SG relative to air, Defaults to 0.75 if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined. If > 0, will change zmethod to 'BNS'
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
        gas_pvt: GasPVT object. If provided, sg/co2/h2s/n2/h2/zmethod/cmethod/tc/pc are extracted from it
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, m, sm3/d). Defaults to False (FIELD)
    """
    if gas_pvt is not None:
        sg = gas_pvt.sg
        co2, h2s, n2, h2 = gas_pvt.co2, gas_pvt.h2s, gas_pvt.n2, gas_pvt.h2
        zmethod, cmethod = gas_pvt.zmethod, gas_pvt.cmethod
        tc, pc = gas_pvt.tc, gas_pvt.pc  # already in oilfield units
    if metric:
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        h = np.asarray(h) * M_TO_FT if not isinstance(h, (int, float)) else h * M_TO_FT
        degf = degc_to_degf(degf)
        r_w = r_w * M_TO_FT
        r_ext = r_ext * M_TO_FT
        if D > 0:
            D = D * D_PER_SM3_TO_D_PER_MSCF  # day/sm3 -> day/Mscf
        if gas_pvt is None:  # tc/pc from gas_pvt are already oilfield units
            if tc > 0:
                tc = tc * 1.8  # K to deg R
            if pc > 0:
                pc = pc * BAR_TO_PSI
    if r_w <= 0:
        raise ValueError("Wellbore radius r_w must be positive")
    if r_ext <= r_w:
        raise ValueError("External radius r_ext must be greater than wellbore radius r_w")

    k, h, pr, pwf = np.asarray(k), np.asarray(h), np.asarray(pr), np.asarray(pwf)
    validate_pe_inputs(p=pr)
    validate_pe_inputs(p=pwf)
    if gas_pvt is None:
        zmethod, cmethod, tc, pc = _prepare_gas_rate_inputs(degf, sg, co2, h2s, n2, h2, zmethod, cmethod, tc, pc)
    direction, delta_mp = _compute_delta_mp(pr, pwf, degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s)

    qg = darcy_gas(delta_mp, k, h, degf, r_w, r_ext, S, D, radial=True)
    result = direction * qg
    if metric:
        return result * MSCF_TO_SM3  # Mscf/d -> sm3/d
    return result

def gas_rate_linear(
    k: npt.ArrayLike,
    pr: npt.ArrayLike,
    pwf: npt.ArrayLike,
    area: npt.ArrayLike,
    length: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    sg: float = 0.75,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    gas_pvt = None,
    metric: bool = False,
) -> np.ndarray:
    """ Returns gas rate for linear flow (mscf/day) using Darcy steady state equation & gas pseudopressure
        k: Permeability (mD)
        pr: Reservoir pressure (psia)
        pwf: BHFP (psia)
        area: Net cross sectional area perpendicular to direction of flow (ft2).
        length: Length over which flow takes place (ft)
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                 defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
               Defaults to 'PMC' if not specified
                sg: Gas SG relative to air, Defaults to 0.75 if not specified
        degf: Reservoir Temperature (deg F).
        co2: Molar fraction of CO2. Defaults to zero if not specified
        h2s: Molar fraction of H2S. Defaults to zero if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if not specified
        h2: Molar fraction of Hydrogen. Defaults to zero if not specified
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
        gas_pvt: GasPVT object. If provided, sg/co2/h2s/n2/h2/zmethod/cmethod/tc/pc are extracted from it
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, m, m2, sm3/d). Defaults to False (FIELD)
    """
    if gas_pvt is not None:
        sg = gas_pvt.sg
        co2, h2s, n2, h2 = gas_pvt.co2, gas_pvt.h2s, gas_pvt.n2, gas_pvt.h2
        zmethod, cmethod = gas_pvt.zmethod, gas_pvt.cmethod
        tc, pc = gas_pvt.tc, gas_pvt.pc  # already in oilfield units
    if metric:
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        area = np.asarray(area) * SQM_TO_SQFT if not isinstance(area, (int, float)) else area * SQM_TO_SQFT
        degf = degc_to_degf(degf)
        length = length * M_TO_FT
        if gas_pvt is None:  # tc/pc from gas_pvt are already oilfield units
            if tc > 0:
                tc = tc * 1.8  # K to deg R
            if pc > 0:
                pc = pc * BAR_TO_PSI
    if length <= 0:
        raise ValueError("Flow length must be positive")

    k, area, pr, pwf = np.asarray(k), np.asarray(area), np.asarray(pr), np.asarray(pwf)
    validate_pe_inputs(p=pr)
    validate_pe_inputs(p=pwf)
    if gas_pvt is None:
        zmethod, cmethod, tc, pc = _prepare_gas_rate_inputs(degf, sg, co2, h2s, n2, h2, zmethod, cmethod, tc, pc)
    direction, delta_mp = _compute_delta_mp(pr, pwf, degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s)

    qg = darcy_gas(delta_mp, k, 1, degf, area, length, 0, 0, radial=False)
    result = direction * qg
    if metric:
        return result * MSCF_TO_SM3  # Mscf/d -> sm3/d
    return result

def darcy_gas(
    delta_mp: npt.ArrayLike,
    k: npt.ArrayLike,
    h: npt.ArrayLike,
    degf: float,
    l1: float,
    l2: float,
    S: float,
    D: float,
    radial: bool,
) -> np.ndarray:
    # Returns mscf/day gas rate. k (mD), h (ft), t (deg F), l1 (r_w or width)/l2 (re or length) (ft), S(Skin), D(Day/mscf)
    validate_pe_inputs(degf=degf)
    tr = degf + degF2R
    if radial:
        a = k * h * delta_mp
        b = _DARCY_GAS_CONST * tr
        c = np.log(l2 / l1) - 0.75 + S
        if D > 1e-9:  # Solve analytically for rate with non-Darcy factor by rearranging into root of a quadratic equation.
            return (np.sqrt(4 * a * b * D + (b * b * c * c)) - (b * c)) / (2 * b * D)
    else:
        a = k * h * l1 * delta_mp
        b = _DARCY_GAS_CONST * tr
        c = l2
    # Else, ignore non-Darcy skin
    return a / (b * c)

def gas_tc_pc(
    sg: float,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    cmethod: str = "PMC",
    tc: float = 0,
    pc: float = 0,
    metric: bool = False,
) -> Tuple:
    """ Returns a Tuple of critical temperature (deg R) and critical pressure (psia).
        For SUT and PMC, this returns an equivalent set of critical parameters for the mixture
        (full gas including inerts). For BNS, this returns critical parameters for the
        *inert-free* hydrocarbon gas only: the supplied inert fractions (co2, h2s, n2, h2)
        are used solely to back out the inert-free hydrocarbon specific gravity from the
        overall mixture sg, and the BNS pseudo-critical correlation is then applied to that
        hydrocarbon sg. The returned Tc/Pc therefore describe the pure hydrocarbon pseudo-
        component; the BNS 5-component PR-EOS handles the inert species' Tc/Pc separately
        via its own per-component internal constants (CO2, H2S, N2, H2).
        sg: Specific gravity of reservoir gas (relative to air)
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        cmethod: 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections,
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
        metric: If True, input/output in Eclipse METRIC units (K, barsa). Defaults to False (FIELD)
    """
    validate_pe_inputs(sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    if metric:
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    if tc * pc > 0:  # Critical properties have both been user specified
        if metric:
            return (tc / 1.8, pc * PSI_TO_BAR)  # deg R -> K, psia -> barsa
        return (tc, pc)
    
    _, cmethod = _h2_method_override(h2, 'DAK', cmethod)
    cmethod = validate_methods(["cmethod"], [cmethod])

    if cmethod.name == "PMC":  # Piper, McCain & Corredor (1999)
        y = np.array([0, h2s, co2, n2])
        alpha = _PMC_ALPHA
        beta = _PMC_BETA
        tci = _PMC_TCI
        pci = _PMC_PCI
        j = alpha[0] + (alpha[4] * sg) + (alpha[5] * sg * sg)  # 2.5
        k = beta[0] + (beta[4] * sg) + (beta[5] * sg * sg)  # 2.6
        jt = j
        kt = k
        jt += sum([(alpha[i] * y[i] * tci[i] / pci[i]) for i in range(1, 4)])
        kt += sum(
            (beta[i] * y[i] * tci[i] / np.sqrt(pci[i])) for i in range(1, 4)
        )
        tpc = kt * kt / jt  # 2.4
        j += sum([alpha[i] * y[i] * tci[i] / pci[i] for i in range(1, 4)])
        k += sum(
            [beta[i] * y[i] * tci[i] / np.sqrt(pci[i]) for i in range(1, 4)]
        )
        ppc = (k * k / j) / j

    elif (cmethod.name == "SUT"):  # Sutton equations with Wichert & Aziz corrections
        hc_frac = 1 - n2 - co2 - h2s
        if hc_frac <= 1e-6:
            raise ValueError("SUT method requires hydrocarbon fraction > 0 (n2 + co2 + h2s must be < 1.0)")
        sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / MW_AIR) / hc_frac  # Eq 3.53
        ppc_hc = _SUT_PPC[0] + _SUT_PPC[1] * sg_hc + _SUT_PPC[2] * sg_hc ** 2  # Eq 3.47b
        tpc_hc = _SUT_TPC[0] + _SUT_TPC[1] * sg_hc + _SUT_TPC[2] * sg_hc ** 2  # Eq 3.47a

        # Wichert & Aziz non-hydrocarbon corrections from monograph
        eps = _SUT_WA[0] * ((co2 + h2s) ** _SUT_WA[1] - (co2 + h2s) ** _SUT_WA[2]) + _SUT_WA[3] * (
            h2s ** _SUT_WA[4] - h2s ** _SUT_WA[5]
        )  # Eq 3.52c
        ppc_star = (
            (1 - n2 - co2 - h2s) * ppc_hc
            + n2 * _IMPUR_PC[0]
            + co2 * _IMPUR_PC[1]
            + h2s * _IMPUR_PC[2]
        )  # Eq 3.54a
        tpc_star = (
            (1 - n2 - co2 - h2s) * tpc_hc
            + n2 * _IMPUR_TC[0]
            + co2 * _IMPUR_TC[1]
            + h2s * _IMPUR_TC[2]
        )  # Eq 3.54b
        tpc = tpc_star - eps  # Eq 3.52a
        ppc = (
            ppc_star * (tpc_star - eps) / (tpc_star + h2s * (1 - h2s) * eps)
        )  # Eq. 3.52b


    elif cmethod.name in ("BUR", "BNS"):
        def tc_ag(x):
            a, b, c = _BNS_TC_AG
            return a * x / (b + x) + c
        
        def tc_gc(x):
            a, b, c = _BNS_TC_GC
            return a * x / (b + x) + c
        
        def pc_fn(x, vc_slope, tc_p):
            vc_on_zc = vc_slope * x + _BNS_VC_INTERCEPT
            return R * tc_p / vc_on_zc
        
        def pseudo_critical(sg_hc, AG = False):
            """
            Calculates pseudo-critical temperature and pressure for the pseudo-hydrocarbon component.
            - Uses custom linear fits for gas condensates (AG = False) and associated gas (AG = True) (Burgoyne, 2025).
            - These relations are derived from property fitting and may differ from literature.
            """
            x = max(0, MW_AIR * sg_hc - _BNS_CH4_MW)
            if AG:
                tpc_hc = tc_ag(x)
                vc_slope = _BNS_VC_SLOPE_AG
            else:
                tpc_hc = tc_gc(x)
                vc_slope = _BNS_VC_SLOPE_GC
            ppc_hc = pc_fn(x, vc_slope, tpc_hc)
            return tpc_hc, ppc_hc   
            
        if co2 + h2s + n2 + h2 < 1.0: # If not 100% Inerts, then calculate hydrocarbon MW
            hydrocarbon_specific_gravity = (sg - (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2) / MW_AIR) / (1 - co2 - h2s - n2 - h2)
        else:
            hydrocarbon_specific_gravity = 0.75 # Use default value if 100% inerts to avoid numerical problems
        hydrocarbon_specific_gravity = np.max([_BNS_SG_METHANE, hydrocarbon_specific_gravity])  # Methane is lower limit
        tpc, ppc = pseudo_critical(hydrocarbon_specific_gravity)

    else:
        raise ValueError("Incorrect cmethod specified")

    if tc > 0:
        tpc = tc
    if pc > 0:
        ppc = pc
    if metric:
        return (tpc / 1.8, ppc * PSI_TO_BAR)  # deg R -> K, psia -> barsa
    return (tpc, ppc)

# EOS parameters for BNS Peng-Robinson model (shared by gas_z and gas_ug)
_BNS_MWS = np.array([44.01, 34.082, 28.014, 2.016, 0])
_BNS_TCS = np.array([547.416, 672.120, 227.160, 47.430, 1])  # H2 Tc has been modified
_BNS_PCS = np.array([1069.51, 1299.97, 492.84, 187.5300, 1])
_BNS_ACF = np.array([0.12253, 0.04909, 0.037, -0.21700, -0.03899])
_BNS_VSHIFT = np.array([-0.27607, -0.22901, -0.21066, -0.36270, -0.19076])
_BNS_OMEGAA = np.array([0.427671, 0.436725, 0.457236, 0.457236, 0.457236])
_BNS_OMEGAB = np.array([0.0696397, 0.0724345, 0.0777961, 0.0777961, 0.0777961])
_BNS_VCVIS = np.array([1.46352, 1.46808, 1.35526, 0.68473, 0.0])  # cuft/lbmol

# BIP precomputed matrices: kij[i,j] = _BIP_CONST[i,j] + _BIP_SLOPE_TC[i,j] / degR
# Component order: [CO2=0, H2S=1, N2=2, H2=3, Gas=4]
# Gas column/row uses tpc_hc at runtime; stored slopes in _BIP_GAS_SLOPES
_BIP_CONST = np.array([
    [ 0.      ,  0.248638, -0.25    , -0.247153, -0.145561],
    [ 0.248638,  0.      , -0.204414,  0.      ,  0.16852 ],
    [-0.25    , -0.204414,  0.      , -0.166253, -0.108   ],
    [-0.247153,  0.      , -0.166253,  0.      , -0.0620119],
    [-0.145561,  0.16852 , -0.108   , -0.0620119,  0.      ]])

_BIP_SLOPE_TC = np.array([
    [  0.        , -75.64467996,  63.51120432,  89.65031832, 0.],
    [-75.64467996,   0.        , 157.55635404,   0.        , 0.],
    [ 63.51120432, 157.55635404,   0.        ,  17.90313836, 0.],
    [ 89.65031832,   0.        ,  17.90313836,   0.        , 0.],
    [  0.        ,   0.        ,   0.        ,   0.        , 0.]])

_BIP_GAS_SLOPES = np.array([0.276572, -0.122378, 0.0605506, 0.0427873])

def _calc_bips_fast(degR, tpc_hc):
    """Compute 5x5 BIP matrix using precomputed constants."""
    slope_tc = _BIP_SLOPE_TC.copy()
    slope_tc[4, :4] = _BIP_GAS_SLOPES * tpc_hc
    slope_tc[:4, 4] = slope_tc[4, :4]
    return _BIP_CONST + slope_tc / degR

def _cardano_cubic(c2, c1, c0, flag=0):
    """Analytic Cardano solver for monic cubic Z^3 + c2*Z^2 + c1*Z + c0 = 0.
    flag=1: max root, flag=-1: min root, flag=0: all real roots."""
    p = (3 * c1 - c2**2) / 3
    q = (2 * c2**3 - 9 * c2 * c1 + 27 * c0) / 27
    root_diagnostic = q**2 / 4 + p**3 / 27

    if root_diagnostic < 0:
        m = 2 * np.sqrt(-p / 3)
        qpm = np.clip(3 * q / p / m, -1.0, 1.0)
        theta1 = np.arccos(qpm) / 3
        roots = np.array([m * np.cos(theta1),
                          m * np.cos(theta1 + 4 * np.pi / 3),
                          m * np.cos(theta1 + 2 * np.pi / 3)])
        Zs = roots - c2 / 3
    else:
        P = (-q / 2 + np.sqrt(root_diagnostic))
        if P >= 0:
            P = P ** (1 / 3)
        else:
            P = -(-P) ** (1 / 3)
        Q = (-q / 2 - np.sqrt(root_diagnostic))
        if Q >= 0:
            Q = Q ** (1 / 3)
        else:
            Q = -(-Q) ** (1 / 3)
        Zs = np.array([P + Q]) - c2 / 3

    if flag == -1:
        return min(Zs)
    if flag == 1:
        return max(Zs)
    return Zs

def _halley_cubic_vec(c2, c1, c0, A=None, B=None, max_iter=50, tol=1e-12):
    """Vectorized Halley solver: solve Z^3+c2*Z^2+c1*Z+c0=0 for max root (vapor Z).
    c2, c1, c0 are 1D arrays of length N. Returns 1D array of Z values.
    When A and B are provided, solves in Z* = Z - B space (per Aaron Zick's
    reformulation) where all physical roots lie in (0, 1], giving more
    robust convergence. Falls back to _cardano_cubic for any bad elements."""

    if A is not None and B is not None:
        # Z* = Z - B reformulation: Z*³ + d2·Z*² + d1·Z* + d0 = 0
        # f*(0) = -2B² < 0, f*(1) = A >= 0, so largest root in (0, 1]
        d2 = 4.0 * B - 1.0
        d1 = A + 2.0 * B * (B - 2.0)
        d0 = -2.0 * B * B

        # Start at Z* = 1 (above largest root since f*(1) = A >= 0)
        Zs = np.ones_like(c2)

        for _ in range(max_iter):
            f = Zs**3 + d2 * Zs**2 + d1 * Zs + d0
            fp = 3.0 * Zs**2 + 2.0 * d2 * Zs + d1
            fpp = 6.0 * Zs + 2.0 * d2
            safe_fp = np.where(np.abs(fp) < 1e-30, 1e-30, fp)
            dZ = f / safe_fp
            denom = safe_fp - 0.5 * dZ * fpp
            denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
            dZ = f / denom
            Zs -= dZ
            if np.max(np.abs(dZ)) < tol:
                break

        # Convert Z* -> Z; fallback to Cardano in Z space for bad elements
        Z = Zs + B
        f = Zs**3 + d2 * Zs**2 + d1 * Zs + d0
        bad = (np.abs(f) > 1e-6) | (Zs <= 0.0)
    else:
        # Legacy path: solve directly in Z space with Cauchy upper bound
        Z = 1.0 + np.maximum(np.abs(c2), np.maximum(np.abs(c1), np.abs(c0)))

        for _ in range(max_iter):
            f = Z**3 + c2 * Z**2 + c1 * Z + c0
            fp = 3.0 * Z**2 + 2.0 * c2 * Z + c1
            fpp = 6.0 * Z + 2.0 * c2
            safe_fp = np.where(np.abs(fp) < 1e-30, 1e-30, fp)
            dZ = f / safe_fp
            denom = safe_fp - 0.5 * dZ * fpp
            denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
            dZ = f / denom
            Z -= dZ
            if np.max(np.abs(dZ)) < tol:
                break

        f = Z**3 + c2 * Z**2 + c1 * Z + c0
        bad = (np.abs(f) > 1e-6) | (Z < 0.0)

    if np.any(bad):
        bad_idx = np.where(bad)[0]
        for idx in bad_idx:
            Z[idx] = _cardano_cubic(c2[idx], c1[idx], c0[idx], flag=1)

    return Z

def gas_z(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: str = "DAK",
    cmethod: str = "PMC",
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns real-gas deviation factor (Z). Returning either single float, or numpy array depending upon
        whether single pressure of list/array or pressures has been specified.
        p: Gas pressure (psia | barsa)
        sg: Gas SG relative to air. Single float only
        degf: Gas Temperature (deg F | deg C)
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                 defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
               Defaults to 'PMC'
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, K). Defaults to False (FIELD)
    """
    p, degf, tc, pc = _metric_to_field_pvt(p, degf, tc, pc, metric)
    validate_pe_inputs(p=p, degf=degf, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)

    tolerance = 1e-6
    p, is_list = convert_to_numpy(p)
    user_tc_pc = tc > 0 and pc > 0

    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)

    tc, pc = gas_tc_pc(sg, co2, h2s, n2, h2, cmethod.name, tc, pc)
    tr = (degf + degF2R) / tc
    pprs = np.array(p/pc)

    # Correlation validity range warnings (non-blocking)
    zname = zmethod.name
    if zname == 'DAK':
        if tr < 1.05 or tr > 3.0:
            warnings.warn(f"DAK Z-factor: Tr={tr:.3f} outside calibration range [1.05, 3.0]", stacklevel=2)
        if np.any(pprs < 0.2) or np.any(pprs > 30):
            warnings.warn(f"DAK Z-factor: Ppr outside calibration range [0.2, 30]", stacklevel=2)
    elif zname == 'HY':
        if tr < 1.15 or tr > 3.0:
            warnings.warn(f"HY Z-factor: Tr={tr:.3f} outside calibration range [1.15, 3.0]", stacklevel=2)
        if np.any(pprs > 24):
            warnings.warn(f"HY Z-factor: Ppr outside calibration range [0, 24]", stacklevel=2)
    elif zname == 'WYW':
        if tr < 1.0 or tr > 3.0:
            warnings.warn(f"WYW Z-factor: Tr={tr:.3f} outside calibration range [1.0, 3.0]", stacklevel=2)
        if np.any(pprs < 0.01) or np.any(pprs > 30):
            warnings.warn(f"WYW Z-factor: Ppr outside calibration range [0.01, 30]", stacklevel=2)

    def zdak(pprs, tr):
        # DAK from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        # Vectorized Newton-Raphson on reduced density rhor
        A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11 = _DAK_A

        R1 = A1 + A2 / tr + A3 / tr**3 + A4 / tr**4 + A5 / tr**5
        R2 = _DAK_LEADING * pprs / tr  # Array
        R3 = A6 + A7 / tr + A8 / tr**2
        R4 = A9 * (A7 / tr + A8 / tr**2)
        R5 = A10 / tr**3

        rhor = _DAK_LEADING * pprs / tr  # Initial guess (array)
        rhor = np.maximum(rhor, 1e-10)

        for _ in range(100):
            r2 = rhor ** 2
            r5 = rhor ** 5
            exp_term = np.exp(-A11 * r2)
            f_val = (R1 * rhor - R2 / rhor + R3 * r2 - R4 * r5 +
                     R5 * r2 * (1 + A11 * r2) * exp_term + 1)
            fp_val = (R1 + R2 / r2 + 2 * R3 * rhor - 5 * R4 * rhor**4 +
                      2 * R5 * rhor * exp_term *
                      ((1 + 2 * A11 * rhor**3) - A11 * r2 * (1 + A11 * r2)))
            fp_val = np.where(np.abs(fp_val) < 1e-30, 1e-30, fp_val)
            # Converge on both function value and step size (matching original scalar logic)
            if np.all(np.abs(f_val) < tolerance):
                break
            step = f_val / fp_val
            rhor_new = rhor - step
            rhor_new = np.maximum(rhor_new, 1e-10)
            if np.all(np.abs(rhor - rhor_new) < tolerance):
                rhor = rhor_new
                break
            rhor = rhor_new

        zout = _DAK_LEADING * pprs / (rhor * tr)
        return process_output(zout, is_list)

    # Hall & Yarborough — Vectorized Newton-Raphson
    def z_hy(pprs, tr):

        tpr_inv = 1/tr  # Reciprocal reduced temperature
        t2 = tpr_inv ** 2
        _a0, _a1, _b0, _b1, _b2, _c0, _c1, _c2, _d0, _d1 = _HY_COEFFS
        a = _a0 * tpr_inv * np.exp(_a1 * (1 - tpr_inv) ** 2)
        b = tpr_inv * (_b0 + _b1 * tpr_inv + _b2 * t2)
        c = tpr_inv * (_c0 + _c1 * tpr_inv + _c2 * t2)
        D = _d0 + _d1 * tpr_inv

        # Initial guess from WYW
        z_init = z_wyw(pprs, tr)
        z_init = np.atleast_1d(z_init).astype(float)
        yi = np.maximum(a * pprs / z_init, 1e-10)
        # Match original scalar: y starts at 0.01, yi starts at WYW guess
        y = np.full_like(yi, 0.01)

        for _ in range(100):
            # Convergence check at top of loop (original while-loop semantics)
            rel_err = np.abs(y - yi) / np.maximum(np.abs(y), 1e-10)
            if np.all(rel_err <= 0.0005):
                break
            # Newton-Raphson step on yi
            yi_safe = np.clip(yi, 1e-10, 0.99)
            omy3 = (1 - yi_safe) ** 3
            omy4 = (1 - yi_safe) ** 4
            f_val = ((yi_safe + yi_safe**2 + yi_safe**3 - yi_safe**4) / omy3) - a * pprs - b * yi_safe**2 + c * yi_safe**D
            df_val = ((1 + 4*yi_safe + 4*yi_safe**2 - 4*yi_safe**3 + yi_safe**4) / omy4) - 2*b*yi_safe + c*D * yi_safe**(D-1)
            df_val = np.where(np.abs(df_val) < 1e-30, 1e-30, df_val)
            y = yi_safe - f_val / df_val
            y = np.maximum(y, 1e-10)
            yi = y

        y = np.maximum(y, 1e-30)
        zout = a * pprs / y
        return process_output(zout, is_list)

    # Wang, Ye & Wu, 2021, 0.2 < Ppr < 30, 1.05 < tpr < 3.0
    # "An accurate correlation for calculating natural gas compressibility factors under a wide range of pressure conditions"
    # https://doi.org/10.1016/j.egyr.2021.11.029
    def z_wyw(pprs, tr):
        a = _WYW_A
        numerators = a[1] + a[2] * (1 + a[3] * tr + a[4] * tr ** 2 + a[5] * tr ** 3 + a[6] * tr ** 4) * pprs + a[7] * pprs ** 2 + a[8] * pprs ** 3 + a[9] * pprs ** 4
        denominators = a[10] + a[11] * (1 + a[12] * tr + a[13] * tr ** 2 + a[14] * tr ** 3 + a[15] * tr ** 4 + a[16] * tr ** 5) * pprs + a[17] * pprs ** 2 + a[18] * pprs ** 3 + a[19] * pprs ** 4 + a[20] * pprs ** 5
        zout = numerators/denominators
        return process_output(zout, is_list)

    mws, tcs, pcs = _BNS_MWS.copy(), _BNS_TCS.copy(), _BNS_PCS.copy()
    ACF, VSHIFT = _BNS_ACF, _BNS_VSHIFT
    OmegaA, OmegaB, VCVIS = _BNS_OMEGAA, _BNS_OMEGAB, _BNS_VCVIS

    # BNS tuned Peng Robinson EOS
    # More information about formulation and applicability can be found here; https://github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor
    def z_bur(psias, degf):
        degR = degf + degF2R

        z = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])

        tcs[-1], pcs[-1] = tc, pc  # Hydrocarbon Tc and Pc from SG using BNS correlation
        trs = degR / tcs

        m = 0.37464 + 1.54226 * ACF - 0.26992 * ACF**2
        alpha = (1 + m * (1 - np.sqrt(trs)))**2

        kij = _calc_bips_fast(degR, tc)

        # Vectorized across all pressures (N = len(psias))
        prs = psias[:, None] / pcs[None, :]            # (N, 5)
        Ai = OmegaA * alpha * prs / trs**2             # (N, 5)
        Bi = OmegaB * prs / trs                        # (N, 5)

        # Mixing rule A: A = sum_ij z_i*z_j*sqrt(Ai_i*Ai_j)*(1-kij)
        sqrt_Ai = np.sqrt(Ai)                          # (N, 5)
        w = z * sqrt_Ai                                # (N, 5)
        onemk = 1.0 - kij                              # (5, 5)
        A = np.sum((w @ onemk) * w, axis=1)            # (N,)

        # Mixing rule B
        B = Bi @ z                                      # (N,)

        # Cubic coefficients: Z^3 + c2*Z^2 + c1*Z + c0 = 0
        c2 = -(1.0 - B)
        c1_coeff = A - 3.0 * B**2 - 2.0 * B
        c0 = -(A * B - B**2 - B**3)

        # Solve all cubics at once - get max (vapor) root
        Z_raw = _halley_cubic_vec(c2, c1_coeff, c0, A=A, B=B)    # (N,)

        # Fugacity-based root selection for sub-critical conditions
        # When 3 real roots exist, the thermodynamically stable phase
        # is the one with the lowest fugacity coefficient (Gibbs criterion).
        # Discriminant of depressed cubic: disc < 0 means 3 real roots
        p_d = (3.0 * c1_coeff - c2**2) / 3.0
        q_d = (2.0 * c2**3 - 9.0 * c2 * c1_coeff + 27.0 * c0) / 27.0
        disc = q_d**2 / 4.0 + p_d**3 / 27.0
        three_roots = disc < -1e-15

        if np.any(three_roots):
            idx = three_roots
            Z_max_s = Z_raw[idx]
            A_s, B_s = A[idx], B[idx]

            # Deflate cubic by Z_max to get quadratic for remaining roots
            b_q = c2[idx] + Z_max_s
            c_q = c1_coeff[idx] + Z_max_s * b_q
            det = np.maximum(b_q**2 - 4.0 * c_q, 0.0)
            sqrt_det = np.sqrt(det)
            Z_min = (-b_q - sqrt_det) / 2.0  # Smallest root

            # PR fugacity coefficient: ln(phi) for root selection
            sqrt2 = np.sqrt(2.0)
            s2p1, s2m1 = 1.0 + sqrt2, sqrt2 - 1.0

            def _ln_phi(Zv):
                return ((Zv - 1.0) - np.log(Zv - B_s)
                        - A_s / (2.0 * sqrt2 * B_s)
                        * np.log((Zv + s2p1 * B_s) / (Zv - s2m1 * B_s)))

            # Only compare where min root is physically valid (Z > B)
            valid = Z_min > B_s
            Z_min_safe = np.where(valid, Z_min, Z_max_s)
            use_min = valid & (_ln_phi(Z_min_safe) < _ln_phi(Z_max_s))
            Z_raw[idx] = np.where(use_min, Z_min, Z_max_s)

        # Volume translation
        vshift = np.sum(z * VSHIFT * Bi, axis=1)       # (N,)
        zout = Z_raw - vshift

        return process_output(zout, is_list)
        
    zfuncs = {"DAK": zdak, "HY": z_hy, "WYW": z_wyw, "BNS": z_bur}

    # Rust acceleration dispatch (batch — single FFI call for all pressures).
    # User tc/pc is forwarded to each Rust batch path; for BNS it overrides only
    # the HC pseudo-component (inert Tc/Pc stay at BNS internal constants).
    tc_arg = float(tc) if user_tc_pc else 0.0
    pc_arg = float(pc) if user_tc_pc else 0.0
    if RUST_AVAILABLE:
        if zmethod.name in ('BNS', 'BUR'):
            zout = np.array(_rust_module.bns_zfactor_batch(
                p.tolist(), float(degf), float(sg),
                float(co2), float(h2s), float(n2), float(h2),
                tc_arg, pc_arg,
            ))
            return process_output(zout, is_list)
        elif zmethod.name == 'DAK':
            if cmethod.name == 'SUT':
                zout = np.array(_rust_module.dak_zfactor_batch(
                    p.tolist(), float(degf), float(sg),
                    float(co2), float(h2s), float(n2),
                    tc_arg, pc_arg,
                ))
            else:
                zout = np.array([_rust_module.dak_zfactor(float(pr_i), float(tr)) for pr_i in pprs])
            return process_output(zout, is_list)
        elif zmethod.name == 'HY':
            if cmethod.name == 'SUT':
                zout = np.array(_rust_module.hy_zfactor_batch(
                    p.tolist(), float(degf), float(sg),
                    float(co2), float(h2s), float(n2),
                    tc_arg, pc_arg,
                ))
            else:
                zout = np.array([_rust_module.hall_yarborough_zfactor(float(pr_i), float(tr)) for pr_i in pprs])
            return process_output(zout, is_list)

    if zmethod.name in ('BNS', 'BUR'):
        return zfuncs["BNS"](p, degf)
    else:
        return zfuncs[zmethod.name](pprs, tr)

def gas_ug(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    zee: float = 0,
    ugz = False,
    metric: bool = False,
) -> np.ndarray:
    """ Returns Gas Viscosity (cP) or Gas Viscosity * Z-Factor
        Uses Lee, Gonzalez & Eakin (1966) Correlation using equations 2.14-2.17 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        Except if 'BNS' Z-Factor method chosen, in which case tuned LBC viscosity model will be used instead

          p: Gas pressure (psia)
          sg: Gas SG relative to air
          degf: Reservoir Temperature (deg F).
          zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                    defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
          cmethod: Method for calculting critical properties
                   'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                   'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                   Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
          zee: Gas Z-Factor. If undefined, will trigger Z-Factor calculation.
          ugz: Boolean flag that if True returns ugZ instead of ug
          metric: If True, input/output in Eclipse METRIC units (barsa, degC). Defaults to False (FIELD)
    """
    if metric:
        p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if isinstance(tc, (int, float)) and tc > 0:
            tc = tc * 1.8  # K to deg R
        if isinstance(pc, (int, float)) and pc > 0:
            pc = pc * BAR_TO_PSI
    zee_provided = not (isinstance(zee, (int, float, bool)) and not zee)
    p, is_list = convert_to_numpy(p)
    zee, _ = convert_to_numpy(zee)
    user_tc_pc = isinstance(tc, (int, float)) and isinstance(pc, (int, float)) and tc > 0 and pc > 0

    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)

    t = degf + degF2R
    m = MW_AIR * sg

    if not zee_provided or not check_2_inputs(zee, p): # Need to calculate Z-Factors if not same length / type as p
         zee = gas_z(p, sg, degf, zmethod, cmethod, co2, h2s, n2, h2, tc, pc)
    
    rho = m * p / (t * zee * R * WDEN)

    # Rust-accelerated viscosity (batch — single FFI call for all pressures).
    # For BNS, user tc/pc is forwarded as HC-only override (inert Tc/Pc stay
    # at BNS internal constants).
    if RUST_AVAILABLE:
        if zmethod.name not in ('BNS', 'BUR'):
            ug_list = np.array(_rust_module.gas_ug_lge_batch(
                p.tolist(), zee.tolist(), sg, degf
            ))
            ug = process_output(ug_list, is_list)
        else:
            tc_arg = float(tc) if user_tc_pc else 0.0
            pc_arg = float(pc) if user_tc_pc else 0.0
            ug_list = np.array(_rust_module.gas_ug_lbc_batch(
                p.tolist(), zee.tolist(), sg, degf, co2, h2s, n2, h2,
                tc_arg, pc_arg,
            ))
            ug = process_output(ug_list, is_list)
        if ugz:
            return process_output(ug * zee, is_list)
        else:
            return process_output(ug, is_list)

    if zmethod.name not in ('BNS', 'BUR'):
        b = _LGE[0] + (_LGE[1] / t) + (_LGE[2] * m)  # 2.16
        c = _LGE[3] - (_LGE[4] * b)  # 2.17
        a = ((_LGE[5] + (_LGE[6] * m)) * np.power(t, 1.5) / (_LGE[7] + (_LGE[8] * m) + t))  # 2.15
        ug = process_output(a * _LGE[9] * np.exp(b * np.power(rho, c)), is_list)  # 2.14
    else:
        # Vectorized LBC viscosity for BNS method
        # From https://wiki.whitson.com/bopvt/visc_correlations/
        mws_lbc = _BNS_MWS.copy()
        tcs_lbc = _BNS_TCS.copy()
        pcs_lbc = _BNS_PCS.copy()
        VCVIS_lbc = _BNS_VCVIS.copy()  # Copy to avoid mutating module-level array

        degR = degf + degF2R
        zi = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])
        if n2 + co2 + h2s + h2 < 1:
            sg_hc = (sg - (co2 * mws_lbc[0] + h2s * mws_lbc[1] + n2 * mws_lbc[2] + h2 * mws_lbc[3]) / MW_AIR) / (1 - co2 - h2s - n2 - h2)
        else:
            sg_hc = 0.75

        sg_hc = max(sg_hc, _BNS_SG_METHANE)
        hc_gas_mw = sg_hc * MW_AIR

        mws_lbc[-1] = hc_gas_mw
        # User-supplied tc/pc override only the hydrocarbon pseudo-component Tc/Pc;
        # gas_tc_pc returns user values directly when both are >0, else BNS correlation.
        tcs_lbc[-1], pcs_lbc[-1] = gas_tc_pc(hc_gas_mw / MW_AIR, cmethod='BNS',
                                             tc=tc if user_tc_pc else 0,
                                             pc=pc if user_tc_pc else 0)
        VCVIS_lbc[-1] = _BNS_VCVIS_HC[0] * (hc_gas_mw - _BNS_CH4_MW) + _BNS_VCVIS_HC[1]

        # Vectorized Stiel-Thodos
        Tr = degR / tcs_lbc
        Tc_K = tcs_lbc * 5.0 / 9.0
        Pc_atm = pcs_lbc / psc
        eta_st = Tc_K**(1.0/6.0) / (mws_lbc**0.5 * Pc_atm**(2.0/3.0))
        ui_low = _ST_LOW[0] * Tr**_ST_LOW[1] / eta_st
        ui_high = _ST_HIGH[0] * np.maximum(_ST_HIGH[1] * Tr + _ST_HIGH[2], 1e-30)**(5.0/8.0) / eta_st
        ui = np.where(Tr <= _ST_TR_THRESH, ui_low, ui_high)

        # Herning-Zippener dilute gas mixture viscosity
        sqrt_mws = np.sqrt(mws_lbc)
        u0_val = np.sum(zi * ui * sqrt_mws) / np.sum(zi * sqrt_mws)

        # LBC mixture parameters
        a_lbc = _LBC_POLY
        rhoc = 1.0 / np.sum(VCVIS_lbc * zi)
        eta_mix = np.abs(np.sum(zi * Tc_K))**(1.0/6.0) / (np.abs(np.sum(zi * mws_lbc))**0.5 * np.abs(np.sum(zi * Pc_atm))**(2.0/3.0))

        # Vectorized over pressures
        zee_arr, _ = convert_to_numpy(zee)
        rhor = p / (zee_arr * R * degR * rhoc)
        lhs = a_lbc[0] + rhor * (a_lbc[1] + rhor * (a_lbc[2] + rhor * (a_lbc[3] + rhor * a_lbc[4])))
        ug = process_output((lhs**4 - _LBC_OFFSET) / eta_mix + u0_val, is_list)
    if ugz:
        return process_output(ug * zee, is_list)
    else:
        return process_output(ug, is_list)

def gas_cg(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    metric: bool = False,
) -> np.ndarray:
    """ Returns gas compressibility (1/psi) using the 'DAK' Dranchuk & Abou-Kassem (1975) Z-Factor &
        Critical property correlation values if not explicitly specified
        If h2 > 0, will use the 'BNS' Z-Factor method
        p: Gas pressure (psia)
        sg: Gas SG relative to air. Defaults to False if undefined
        pwf: BHFP (psia)
        degf: Gas Temperature (deg F)
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
        zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                    defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculating Tc and Pc
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                 Defaults to 'PMC
          metric: If True, input/output in Eclipse METRIC units (barsa, degC, 1/barsa). Defaults to False (FIELD)
    """
    p, degf, tc, pc = _metric_to_field_pvt(p, degf, tc, pc, metric)
    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)

    p, is_list = convert_to_numpy(p)
    tc, pc = gas_tc_pc(sg=sg, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc, cmethod=cmethod)

    pr = p / pc
    degR = (degf + degF2R)
    tr = degR / tc
    dp = np.maximum(p * 1e-4, 0.01)  # Relative step, floor at 0.01 psi
    p_both = np.concatenate([p, p + dp])
    zee_both = gas_z(p=p_both, sg=sg, degf=degf, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2=h2, tc=tc, pc=pc)
    n = len(p)
    zee1 = zee_both[:n]
    zee2 = zee_both[n:]

    vol1 = zee1*R*degR/p
    vol2 = zee2*R*degR/(p+dp)

    result = process_output((vol1 - vol2)/((vol1 + vol2)/2) / dp, is_list) # 1/psi
    if metric:
        return result * INVPSI_TO_INVBAR  # 1/psi -> 1/barsa
    return result

def gas_bg(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns Bg (gas formation volume factor) for natural gas (rcf/scf)
        p: Gas pressure (psia)
        sg: Gas SG relative to air
        degf: Reservoir Temperature (deg F)
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                 defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                 Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Nitrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R | K). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
          pc: Critical gas pressure (psia | barsa). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
          metric: If True, input/output in Eclipse METRIC units (barsa, degC). Defaults to False (FIELD)
    """
    p, degf, tc, pc = _metric_to_field_pvt(p, degf, tc, pc, metric)
    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)
    p, is_list = convert_to_numpy(p)

    zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)
    degR = (degf + degF2R)
    return process_output(zee * degR / (p * (tsc + degF2R) / psc), is_list)

def gas_sg(hc_mw: float, co2: float, h2s: float, n2: float, h2: float)  -> float:
    """ Returns sg of gas mixture """
    return (hc_mw * (1 - co2 - h2s-  n2 - h2) + (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2)) / MW_AIR

def gas_den(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns gas density for natural gas (lb/cuft)
          p: Gas pressure (psia)
          sg: Gas SG relative to air
          degf: Reservoir Temperature (deg F)
          zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                   defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
          cmethod: Method for calculting critical properties
                   'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                   'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                   Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R | K). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
          pc: Critical gas pressure (psia | barsa). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
          metric: If True, input/output in Eclipse METRIC units (barsa, degC, kg/m3). Defaults to False (FIELD)
    """
    p, degf, tc, pc = _metric_to_field_pvt(p, degf, tc, pc, metric)
    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)
    p, is_list = convert_to_numpy(p)

    zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)

    m = sg * MW_AIR

    if co2 == 1:
        m = 44.01
    if h2s == 1:
        m = 34.082
    if n2 == 1:
        m = 28.014
    if h2 == 1:
        m = 2.016
    degR = degf + degF2R
    rhog = p * m / (zee * R * degR)

    result = process_output(rhog, is_list)
    if metric:
        return result * LBCUFT_TO_KGM3  # lb/cuft -> kg/m3
    return result

def gas_ponz2p(
    poverz: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    rtol: float = 1e-7,
    metric: bool = False,
) -> np.ndarray:
    """ Returns pressure corresponding to a P/Z value for natural gas (psia)
        Calculated through iterative solution method
        poverz: Gas pressure / Z-Factor (psia)
        sg: Gas SG relative to air
        degf: Reservoir Temperature (deg F).
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                 defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                 Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
          rtol: Relative solution tolerance. Will iterate until abs[(p - poverz * Z)/p] < rtol
          metric: If True, input/output in Eclipse METRIC units (barsa, degC). Defaults to False (FIELD)
    """
    if metric:
        poverz = np.asarray(poverz) * BAR_TO_PSI if not isinstance(poverz, (int, float)) else poverz * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI

    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)

    poverz, is_list = convert_to_numpy(poverz)

    def _ponz2p_err_msg(target):
        return (
            f"gas_ponz2p: no single-phase solution exists for P/Z={target:.4f} at this "
            f"temperature and composition. Target may fall in the two-phase region where "
            f"P/Z is discontinuous."
        )

    # Try Rust-accelerated batch solver for supported method combinations.
    # User tc/pc is honored by the Rust path: mixture override for DAK/HY+SUT,
    # HC-only override for BNS.
    zname = zmethod.name
    cname = cmethod.name
    if cname in ('SUT', 'BNS'):
        try:
            p_list = _rust_module.gas_ponz2p_rust(
                [float(v) for v in poverz], degf, sg,
                zname, cname,
                co2, h2s, n2, h2, tc, pc, rtol,
            )
            result = process_output(np.array(p_list), is_list)
            if metric:
                return result * PSI_TO_BAR
            return result
        except (ImportError, AttributeError):
            pass
        except ValueError as e:
            if "failed to converge" in str(e) or "root not bracketed" in str(e):
                raise ValueError(_ponz2p_err_msg(poverz[0])) from e
            raise

    # Python fallback: scalar bisection
    def PonZ2P_err(args, p):
        ponz, sg, degf, zmethod, cmethod, tc, pc, co2, h2s, n2, h2 = args
        zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)
        return (p - (ponz * zee)) / p

    p = []
    for ponz in poverz:
        args = (ponz, sg, degf, zmethod, cmethod, tc, pc, co2, h2s, n2, h2)
        try:
            p.append(bisect_solve(args, PonZ2P_err, ponz * 0.1, ponz * 5, rtol))
        except (ValueError, RuntimeError) as e:
            raise ValueError(_ponz2p_err_msg(ponz)) from e

    result = process_output(p, is_list)
    if metric:
        return result * PSI_TO_BAR  # psia -> barsa
    return result

def gas_grad2sg(
    grad: float,
    p: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    rtol: float = 1e-7,
    metric: bool = False,
) -> float:
    """ Returns insitu gas specific gravity consistent with observed gas gradient. Solution iteratively calculated via bisection
        Calculated through iterative solution method. Bisection bounds span pure H2 (SG ~0.070) to 3.0 to
        accommodate H2-blend and CO2-rich compositions; results outside this range will fail.

        grad: Observed gas gradient (psi/ft)
        p: Pressure at observation (psia)
        degf: Reservoir Temperature (deg F). Defaults to False if undefined
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                 defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                 Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
          rtol: Relative solution tolerance. Will iterate until abs[(grad - calculation)/grad] < rtol
          metric: If True, input/output in Eclipse METRIC units (bar/m, barsa, degC). Defaults to False (FIELD)
    """
    if metric:
        grad = grad * BARM_TO_PSIFT  # bar/m -> psi/ft
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI

    validate_pe_inputs(p=p, degf=degf, co2=co2, h2s=h2s, n2=n2, h2=h2)
    if grad <= 0:
        raise ValueError("Gas gradient must be positive")

    degR = degf + degF2R

    def grad_err(args, sg):
        grad, p, zmethod, cmethod, tc, pc, co2, h2s, n2, h2 = args
        m = sg * MW_AIR
        zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)
        grad_calc = p * m / (zee * R * degR) / 144
        error = (grad - grad_calc) / grad
        return error

    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)

    args = (grad, p, zmethod, cmethod, tc, pc, co2, h2s, n2, h2)
    return bisect_solve(args, grad_err, _GRAD2SG_SG_LO, _GRAD2SG_SG_HI, rtol)

def gas_dmp(
    p1: float,
    p2: float,
    degf: float,
    sg: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    metric: bool = False,
) -> float:
    """ Numerical integration of real-gas pseudopressure between two pressures
        Returns integral over range between p1 to p2 (psi**2/cP)
        p1: Starting (lower) pressure (psia)
        p2: Ending (upper) pressure (psia)
        t: Gas Temperature (deg F)
        sg: Specific gravity of  gas (relative to air)
        zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'BNS' Tuned 5 component Peng Robinson EOS model, Burgoyne, Nielsen and Stanko (2025), SPE-229932-MS
                   defaults to 'DAK' if not specified, or to 'BNS' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
                 Defaults to 'PMC'
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        tc: Critical gas temperature (deg R | K). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
        pc: Critical gas pressure (psia | barsa). Calculates using cmethod if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, bar^2/cP). Defaults to False (FIELD)
    """
    if metric:
        p1 = p1 * BAR_TO_PSI
        p2 = p2 * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    if p1 < 0:
        raise ValueError(f"Pressure p1 must be non-negative, got {p1}")
    if p2 < 0:
        raise ValueError(f"Pressure p2 must be non-negative, got {p2}")
    if p1 == p2:
        return 0

    zmethod, cmethod = _resolve_methods(zmethod, cmethod, h2=h2)

    # Rust-accelerated pseudopressure (entire integration in Rust — no Python round-trips).
    # User tc/pc is honored by the Rust path: mixture override for DAK/HY+SUT,
    # HC-only override for BNS.
    if RUST_AVAILABLE:
        zname = zmethod.name
        cname = cmethod.name
        if cname in ('SUT', 'BNS'):
            try:
                result = _rust_module.gas_dmp_rust(
                    float(p1), float(p2), degf, sg,
                    zname, cname,
                    co2, h2s, n2, h2, tc, pc,
                )
                if metric:
                    return result * PSI2CP_TO_BAR2CP
                return result
            except ValueError:
                pass  # Fall through to Python for unsupported method combos

    def _gl_integrate(lo, hi, nodes, weights):
        """Batch Gauss-Legendre integration of 2p/(mu*Z) over [lo, hi]."""
        p_mid = (lo + hi) * 0.5
        p_half = (hi - lo) * 0.5
        p_eval = p_mid + p_half * nodes
        zee = gas_z(p=p_eval, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod,
                     co2=co2, h2s=h2s, n2=n2, h2=h2, tc=tc, pc=pc)
        mugz = gas_ug(p_eval, sg, degf, zmethod, cmethod, co2, h2s, n2, h2, tc, pc, zee, ugz=True)
        return p_half * np.sum(weights * 2.0 * p_eval / mugz)

    # Two-tier integration: compute with n=7 and n=10, compare for convergence
    result_7 = _gl_integrate(p1, p2, _GL7_NODES, _GL7_WEIGHTS)
    result_10 = _gl_integrate(p1, p2, _GL10_NODES, _GL10_WEIGHTS)

    if abs(result_10) < 1e-30 or abs(result_10 - result_7) / abs(result_10) < 1e-5:
        result = result_10
    else:
        # If not converged, split into two subintervals and integrate with n=10 each (batch)
        p_mid = (p1 + p2) * 0.5
        p_half_lo = (p_mid - p1) * 0.5
        p_half_hi = (p2 - p_mid) * 0.5
        p_center_lo = (p1 + p_mid) * 0.5
        p_center_hi = (p_mid + p2) * 0.5

        # Build all evaluation points for both subintervals in a single array
        p_eval = np.concatenate([p_center_lo + p_half_lo * _GL10_NODES,
                                 p_center_hi + p_half_hi * _GL10_NODES])

        zee = gas_z(p=p_eval, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod,
                     co2=co2, h2s=h2s, n2=n2, h2=h2, tc=tc, pc=pc)
        mugz = gas_ug(p_eval, sg, degf, zmethod, cmethod, co2, h2s, n2, h2, tc, pc, zee, ugz=True)
        integrand = 2.0 * p_eval / mugz

        n = len(_GL10_NODES)
        result = (p_half_lo * np.sum(_GL10_WEIGHTS * integrand[:n]) +
                  p_half_hi * np.sum(_GL10_WEIGHTS * integrand[n:]))

    if metric:
        return result * PSI2CP_TO_BAR2CP  # psi^2/cP -> bar^2/cP
    return result

def gas_fws_sg(sg_g: float, cgr: float, api_st: float, metric: bool = False) -> float:
    """
     Estimates FWS specific gravity of gas-condensate from separator gas SG, CGR and API
     Uses Standing correlation to estimate condensate MW from API.
     Returns SG of FWS gas

     sg_g: Specific gravity of weighted average surface gas (relative to air)
     api_st: Density of stock tank liquid (API)
     cgr: Condensate gas ratio (stb/mmscf | sm3/sm3 if metric=True)
     metric: If True, cgr is in sm3/sm3. Defaults to False (FIELD: stb/mmscf)
     """
    if metric:
        cgr = cgr / STB_PER_MMSCF_TO_SM3_PER_SM3  # sm3/sm3 -> stb/mmscf
    # 1 mmscf separator gas basis with 379.482 scf/lb-mole
    cond_vol = cgr * CUFTperBBL  # cuft/mmscf surface gas
    cond_sg = 141.5 / (api_st + 131.5)
    cond_mass = cond_vol * (cond_sg * WDEN)  # lb
    surface_gas_moles = 1e6 / scf_per_mol  # lb-moles
    surface_gas_mass = sg_g * MW_AIR * surface_gas_moles  # lb
    cond_mw = _STANDING_COND_MW[0] - _STANDING_COND_MW[1] * api_st  # lb/lb-moles (Standing correlation)
    cond_moles = cond_mass / cond_mw  # lb-moles
    fws_gas_mass = cond_mass + surface_gas_mass  # lb
    fws_gas_moles = cond_moles + surface_gas_moles  # lb-moles
    fws_gas_mw = fws_gas_mass / fws_gas_moles  # lb/lb-moles
    return fws_gas_mw / MW_AIR


def gas_hvf_beta(k: float, method: str = 'FK', phi: float = 0.0, *, metric: bool = False) -> float:
    """ Returns the Forchheimer high-velocity-flow (inertial) coefficient β.

        Three correlations for β as a function of rock properties:
        'FK' -- Log-log fit of the Firoozabadi & Katz (1979) consolidated-rock
                β(k) chart: β[1/ft] = 2.172e10 * k^-1.201, k in md. Default.
                Source data: Firoozabadi, A. & Katz, D.L. "An Analysis of
                High-Velocity Gas Flow Through Porous Media," JPT Feb 1979,
                pp.211-216 (SPE-6827).
        'JONES' -- Jones, S.C. (1987) "Using the Inertial Coefficient b to
                   Characterize Heterogeneity in Reservoir Rock," SPE-16949:
                   β[1/ft] = 6.15e10 * k^-1.55, k in md.
        'TCK' -- Tek, Coats, Katz (1962) "The Effect of Turbulence on Flow of
                 Natural Gas through Porous Reservoirs," JPT July 1962,
                 pp.799-806: β[1/ft] = 1.88e10 / (k^1.47 * φ^0.53).
                 Requires phi > 0.

        k: Permeability (md). To evaluate β at a damaged-zone permeability,
           pass k' = k * krg (absolute permeability times gas krg at critical
           oil saturation).
        method: 'FK' | 'JONES' | 'TCK'. Defaults to 'FK'.
        phi: Porosity fraction, required for 'TCK'.
        metric: If True, returns β in 1/m; else 1/ft. Defaults to False.
    """
    if k <= 0:
        raise ValueError(f"Permeability must be positive, got k={k}")
    m = method.upper()
    if m not in _BETA_METHODS:
        raise ValueError(f"beta method must be one of {_BETA_METHODS}, got '{method}'")
    if m == 'FK':
        beta_ft = _FK_A * k ** (-_FK_N)
    elif m == 'JONES':
        beta_ft = _JONES_A * k ** (-_JONES_N)
    else:  # TCK
        if phi <= 0 or phi >= 1:
            raise ValueError(f"TCK requires 0 < phi < 1, got phi={phi}")
        beta_ft = _TCK_A / (k ** _TCK_NK * phi ** _TCK_NPHI)
    if metric:
        return beta_ft * M_TO_FT  # 1/ft * ft/m = 1/m
    return beta_ft


def gas_non_darcy_skin(qg: float, k: float, h_perf: float, rw: float,
                       mug: float, sg: float,
                       krg: float = 1.0, beta_method: str = 'FK',
                       phi: float = 0.0, *, metric: bool = False) -> dict:
    """ Returns rate-dependent (high-velocity-flow) skin for a gas well.

        Computes β via the selected correlation, the non-Darcy coefficient
        D (Jones 1987 SPE-16949; see also Odeh, Moreland & Schueler 1975
        JPT Dec 1975, pp.1501-1504), and S_hvf = D * qg.

        qg: Gas rate (MSCF/D | sm3/D).
        k: Absolute permeability (md).
        h_perf: Perforated interval thickness (ft | m).
        rw: Wellbore radius (ft | m).
        mug: Gas viscosity at reservoir conditions (cP). Obtain from gas_ug
             or GasPVT.viscosity().
        sg: Gas specific gravity relative to air.
        krg: Gas relative permeability at critical oil saturation. If < 1.0,
             β is evaluated at the damaged-zone permeability k' = k * krg,
             giving a more pessimistic S_hvf estimate. Defaults to 1.0
             (undamaged).
        beta_method: β correlation -- see gas_hvf_beta. Defaults to 'FK'.
        phi: Porosity fraction. Required for beta_method='TCK'.
        metric: If True, inputs in (sm3/D, m, m) and D returned in day/sm3;
                else (MSCF/D, ft, ft) and D in day/MSCF. Defaults to False.

        Returns dict:
          'beta':  Forchheimer coefficient (1/ft or 1/m)
          'D':     non-Darcy coefficient (day/MSCF or day/sm3)
          'S_hvf': rate-dependent skin (dimensionless)
    """
    if qg <= 0:
        raise ValueError(f"Gas rate qg must be positive, got {qg}")
    if k <= 0 or h_perf <= 0 or rw <= 0 or mug <= 0 or sg <= 0:
        raise ValueError("k, h_perf, rw, mug, sg must all be positive")
    if not (0 < krg <= 1.0):
        raise ValueError(f"krg must be in (0, 1.0], got {krg}")
    validate_pe_inputs(sg=sg)

    # Internally work in oilfield units: MSCF/D, ft, ft, md, cp
    if metric:
        qg_mscf = qg * (1.0 / MSCF_TO_SM3)  # sm3/D -> MSCF/D
        h_ft = h_perf * M_TO_FT
        rw_ft = rw * M_TO_FT
    else:
        qg_mscf = qg
        h_ft = h_perf
        rw_ft = rw

    k_eff = k * krg  # damaged-zone permeability for β
    beta_ft = gas_hvf_beta(k_eff, method=beta_method, phi=phi, metric=False)

    # D in day/MSCF uses absolute permeability k; only β sees the optional
    # damaged-zone k' = k*krg.
    D_mscf = _D_COATS_KATZ * beta_ft * sg * k / (mug * h_ft * rw_ft)
    S_hvf = D_mscf * qg_mscf

    if metric:
        beta_out = beta_ft * M_TO_FT  # 1/ft -> 1/m
        D_out = D_mscf / D_PER_SM3_TO_D_PER_MSCF  # day/MSCF -> day/sm3
    else:
        beta_out = beta_ft
        D_out = D_mscf

    return {'beta': beta_out, 'D': D_out, 'S_hvf': S_hvf}


def gas_partial_penetration_skin(htot: float, htop: float, hbot: float,
                                 rw: float, kh_kv: float = 10.0) -> float:
    """ Returns partial-penetration pseudoskin S_p using the analytical
        series solution of Streltsova-Adams, T.D. (1979) "Pressure Drawdown
        in a Well with Limited Flow Entry," SPE J. Nov 1979, pp.1469-1476
        (SPE-7486). Bessel K₀ evaluated via scipy.special.k0; series summed
        until three consecutive relative increments fall below 1e-4 (with
        the symmetry optimisation that only odd n contribute when the
        perforated interval is centrally located).

        All thickness and radius inputs must use a consistent length unit
        (ft or m; only ratios enter the formula).

        htot: Total formation thickness (no-flow-to-no-flow).
        htop: Distance from formation top to top of perforated interval.
        hbot: Distance from formation top to bottom of perforated interval.
        rw: Wellbore radius.
        kh_kv: Horizontal-to-vertical permeability anisotropy ratio (k_h/k_v).
               Defaults to 10.0. Set to 0 to treat vertical permeability as
               negligible (returns S_p = 0, i.e. reservoir behaves as the
               product k*h_perf and no partial-penetration penalty exists).
    """
    from scipy.special import k0 as _k0
    if htot <= 0 or rw <= 0:
        raise ValueError("htot and rw must be positive")
    if not (0.0 <= htop < hbot <= htot):
        raise ValueError(
            f"Must satisfy 0 <= htop < hbot <= htot; got htop={htop}, hbot={hbot}, htot={htot}")
    if kh_kv < 0:
        raise ValueError(f"kh_kv must be >= 0, got {kh_kv}")
    if kh_kv == 0:
        return 0.0  # No vertical permeability => no partial-penetration penalty

    htD = htop / htot
    hbD = hbot / htot
    rD = math.sqrt(1.0 / kh_kv) * rw / htot

    # Vectorised series sum. Drop terms where n·π·rD exceeds the K_0 underflow
    # threshold -- those contributions are negligible anyway.
    n = np.arange(1, _SP_MAX_N + 1)
    x = n * math.pi * rD
    finite = x < 700.0
    n_f = n[finite]
    x_f = x[finite]
    diff = np.sin(n_f * math.pi * hbD) - np.sin(n_f * math.pi * htD)
    terms = diff * diff * _k0(x_f) / (n_f * n_f)
    total = float(terms.sum())

    # Convergence diagnostic: last 5% of terms should contribute negligibly
    tail_start = max(1, int(len(terms) * 0.95))
    tail = float(terms[tail_start:].sum())
    if total != 0 and abs(tail / total) > _SP_CONV_TOL:
        warnings.warn(
            f"gas_partial_penetration_skin: series may not be fully converged at "
            f"N={_SP_MAX_N} (tail/total={tail/total:.2e}). Result may be accurate "
            f"to only ~{abs(tail/total)*100:.1f}%.",
            stacklevel=2)

    return 2.0 / (math.pi ** 2 * (hbD - htD) ** 2) * total


class GasPVT:
    """ Gas PVT wrapper that stores composition and method choices, pre-computes critical properties.
        Exposes z(), viscosity(), density(), bg() methods delegating to gas_z, gas_ug, gas_den, gas_bg.

        sg: Gas SG relative to air. Defaults to 0.75
        co2: Molar fraction of CO2. Defaults to zero
        h2s: Molar fraction of H2S. Defaults to zero
        n2: Molar fraction of Nitrogen. Defaults to zero
        h2: Molar fraction of Hydrogen. Defaults to zero
        zmethod: Method for calculating Z-Factor. Defaults to 'DAK', or 'BNS' if h2 > 0
        cmethod: Method for calculating critical properties. Defaults to 'PMC'
        tc: Critical gas temperature (deg R, or K if metric=True). Defaults to 0 = compute from cmethod.
            For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants).
            Both tc and pc must be > 0 to take effect.
        pc: Critical gas pressure (psia, or barsa if metric=True). Defaults to 0 = compute from cmethod.
            For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants).
            Both tc and pc must be > 0 to take effect.
        metric: If True, methods accept/return Eclipse METRIC units (barsa, degC, kg/m3). Defaults to False

        BNS coupling: if either zmethod or cmethod is BNS, both are forced to BNS for
        thermodynamic consistency (a UserWarning is emitted if this overrules a
        user-chosen non-BNS method). h2 > 0 auto-selects BNS silently.
    """
    def __init__(self, sg=0.75, co2=0, h2s=0, n2=0, h2=0,
                 zmethod='DAK', cmethod='PMC', tc=0, pc=0, metric=False):
        self.sg = sg
        self.co2 = co2
        self.h2s = h2s
        self.n2 = n2
        self.h2 = h2
        self.metric = metric
        # Convert user-supplied metric tc/pc to oilfield units for internal storage
        if metric:
            if tc > 0:
                tc = tc * 1.8  # K -> deg R
            if pc > 0:
                pc = pc * BAR_TO_PSI
        self._user_tc_pc = tc > 0 and pc > 0
        self.zmethod, self.cmethod = _resolve_methods(zmethod, cmethod, h2=h2)
        # self.tc/self.pc reflect the effective critical properties used by methods:
        # user-supplied values when both tc/pc > 0, otherwise cmethod correlation output.
        # For SUT/PMC these are mixture pseudo-critical values; for BNS they are the
        # inert-free hydrocarbon pseudo-critical values.
        self.tc, self.pc = gas_tc_pc(sg, co2, h2s, n2, h2, self.cmethod.name, tc, pc)

    def _convert_inputs(self, p, degf):
        """Convert metric inputs to oilfield for internal calculations."""
        if self.metric:
            p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
            degf = degc_to_degf(degf)
        return p, degf

    def z(self, p, degf):
        """ Returns Z-factor at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        tc = self.tc if self._user_tc_pc else 0
        pc = self.pc if self._user_tc_pc else 0
        return gas_z(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                     cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                     n2=self.n2, h2=self.h2, tc=tc, pc=pc)

    def viscosity(self, p, degf):
        """ Returns gas viscosity (cP) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        tc = self.tc if self._user_tc_pc else 0
        pc = self.pc if self._user_tc_pc else 0
        return gas_ug(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                      cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                      n2=self.n2, h2=self.h2, tc=tc, pc=pc)

    def density(self, p, degf):
        """ Returns gas density (lb/cuft | kg/m3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        tc = self.tc if self._user_tc_pc else 0
        pc = self.pc if self._user_tc_pc else 0
        result = gas_den(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                         cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                         n2=self.n2, h2=self.h2, tc=tc, pc=pc)
        if self.metric:
            return result * LBCUFT_TO_KGM3
        return result

    def bg(self, p, degf):
        """ Returns gas FVF Bg (rcf/scf | rm3/sm3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        tc = self.tc if self._user_tc_pc else 0
        pc = self.pc if self._user_tc_pc else 0
        return gas_bg(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                      cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                      n2=self.n2, h2=self.h2, tc=tc, pc=pc)

    def non_darcy_skin(self, qg, p, degf, k, h_perf, rw,
                       krg=1.0, beta_method='FK', phi=0.0):
        """ Rate-dependent (HVF) skin for this gas. μ_g is computed internally
            from stored PVT state at (p, degf). See gas_non_darcy_skin().

            qg: Gas rate (MSCF/D | sm3/D)
            p: Reservoir pressure (psia | barsa)
            degf: Reservoir temperature (deg F | deg C)
            k, h_perf, rw: Permeability (md), perforated thickness (ft | m),
                           wellbore radius (ft | m).
            krg, beta_method, phi: see gas_non_darcy_skin().

            Returns dict with 'beta', 'D', 'S_hvf'.
        """
        mug = self.viscosity(p, degf)
        return gas_non_darcy_skin(qg=qg, k=k, h_perf=h_perf, rw=rw, mug=mug,
                                  sg=self.sg, krg=krg, beta_method=beta_method,
                                  phi=phi, metric=self.metric)

    def partial_penetration_skin(self, htot, htop, hbot, rw, kh_kv=10.0):
        """ Partial-penetration skin (Streltsova-Adams 1979). No PVT state is
            used; delegated for API consistency. See
            gas_partial_penetration_skin().
        """
        return gas_partial_penetration_skin(htot=htot, htop=htop, hbot=hbot,
                                            rw=rw, kh_kv=kh_kv)


def gas_water_content(p: float, degf: float, salinity: float = 0, metric: bool = False) -> float:
    """ Returns saturated volume of water vapor in natural gas (stb/mmscf | sm3/sm3 if metric)
        From 'PVT and Phase Behaviour Of Petroleum Reservoir Fluids' by Ali Danesh
        degf: Water Temperature (deg F | deg C)
        p: Water pressure (psia | barsa)
        salinity: Water salinity (wt% NaCl). Defaults to 0 (freshwater)
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
    validate_pe_inputs(p=p, degf=degf)
    t = degf
    content = (
        (
            _DANESH_WC[0]
            * (
                np.exp(
                    _DANESH_WC[1]
                    + (_DANESH_WC[2] / (t + degF2R))
                    + (_DANESH_WC[3] * np.log(t + degF2R))
                    + (_DANESH_WC[4] * ((t + degF2R) * (t + degF2R)))
                )
            )
            / (p)
            + (np.power(10, ((_DANESH_WC[5] / (t + degF2R)) + _DANESH_WC[6])))
        )
        * (1 - (_DANESH_WC[7] * salinity) - (_DANESH_WC[8] * (salinity * salinity)))
        / _DANESH_WC[9]
        / _DANESH_WC[10]
    )
    if metric:
        return content * STB_PER_MMSCF_TO_SM3_PER_SM3  # stb/MMscf -> sm3/sm3
    return content


# =============================================================================
# Gas Hydrate Formation Prediction
# =============================================================================

_PSI_TO_KPA = 6.894757  # psia to kPa

# Hammerschmidt (1934) inhibitor constants: (K, MW, w_max as fraction)
_HAMMERSCHMIDT = {
    inhibitor.MEOH: (2335, 32.04, 0.25),
    inhibitor.MEG:  (2700, 62.07, 0.70),
    inhibitor.DEG:  (2700, 106.12, 0.70),
    inhibitor.TEG:  (5400, 150.17, 0.50),
    inhibitor.ETOH: (1297, 46.07, 0.30),
}

# Østergaard et al. (2005) coefficients: delta_T(degC) = C1*w + C2*w^2 + C3*w^3, w in wt% (0-100)
_OSTERGAARD = {
    inhibitor.MEOH: (0.4411, -0.0033, 6.476e-5),
    inhibitor.MEG:  (0.2533, -0.0009, 2.222e-5),
    inhibitor.DEG:  (0.1833, -0.0004, 6.667e-6),
    inhibitor.TEG:  (0.1333, -0.0002, 3.333e-6),
    inhibitor.ETOH: (0.3750, -0.0020, 3.500e-5),
}


# Inhibitor physical density at 20 degC (g/cm3) — for volumetric injection rate
_INHIBITOR_DENSITY = {
    inhibitor.MEOH: 0.791,
    inhibitor.MEG:  1.110,
    inhibitor.DEG:  1.117,
    inhibitor.TEG:  1.125,
    inhibitor.ETOH: 0.789,
}

# Conversion: g/cm3 to lb/gal
_GCM3_TO_LB_PER_GAL = 8.34540445

# Water mass at standard conditions: lb per stb
_WATER_LB_PER_STB = 350.2

# Maximum valid wt% for each inhibitor (aqueous phase concentration)
_MAX_WT_PCT = {
    inhibitor.MEOH: 25.0,
    inhibitor.MEG:  70.0,
    inhibitor.DEG:  70.0,
    inhibitor.TEG:  50.0,
    inhibitor.ETOH: 30.0,
}


@dataclass
class HydrateResult:
    """Result of gas hydrate formation prediction.

    Hydrate assessment (HFT, HFP, subcooling, inhibitor) is evaluated at the
    operating point (p, degf). Water balance is evaluated between reservoir
    conditions (p_res, degf_res) and the operating point.

    Attributes — Hydrate Assessment
    --------------------------------
    hft : float
        Hydrate formation temperature at operating pressure (degF | degC).
    hfp : float
        Hydrate formation pressure at operating temperature (psia | barsa).
    subcooling : float
        HFT - T_operating (degF | degC delta). Positive means inside hydrate window.
    in_hydrate_window : bool
        True if operating temperature is below the hydrate formation temperature.
    inhibited_hft : float
        HFT after inhibitor depression (degF | degC), or NaN if no inhibitor.
    inhibitor_depression : float
        Temperature depression from inhibitor (degF | degC delta), or 0.
    required_inhibitor_wt_pct : float
        Wt% inhibitor in aqueous phase (water + inhibitor) needed to bring HFT
        below operating temperature, or 0. Capped at the physical maximum for
        the selected inhibitor type.
    max_inhibitor_wt_pct : float
        Maximum valid wt% for the selected inhibitor type (MEOH: 25%, MEG: 70%,
        DEG: 70%, TEG: 50%, ETOH: 30%). 0 if no inhibitor specified.
    inhibitor_underdosed : bool
        True if required_inhibitor_wt_pct exceeds max_inhibitor_wt_pct,
        meaning this inhibitor type cannot provide sufficient depression
        even at its physical maximum concentration. Does NOT indicate
        whether the applied inhibitor_wt_pct is sufficient — compare
        inhibited_hft to operating temperature to check applied-dose
        protection.

    Attributes — Water Balance
    --------------------------
    The gas leaves the reservoir saturated with vaporized water at reservoir
    P,T. At the operating point (lower P,T), the gas can hold less water vapor,
    so the excess condenses as liquid. Free water is any additional liquid water
    entrained in the gas stream from the reservoir (user-specified).

    water_vaporized_res : float
        Equilibrium vaporized water in the gas at reservoir P,T
        (stb/MMscf | sm3/sm3). When p_res/degf_res not provided, equals
        water_vaporized_op (no temperature/pressure change, no condensation).
    water_vaporized_op : float
        Equilibrium vaporized water in the gas at operating P,T
        (stb/MMscf | sm3/sm3).
    water_condensed : float
        Water condensed from vapor between reservoir and operating conditions
        = max(water_vaporized_res - water_vaporized_op, 0) (stb/MMscf | sm3/sm3).
    free_water : float
        Free liquid water influx from reservoir (= additional_water input)
        (stb/MMscf | sm3/sm3).
    total_liquid_water : float
        Total liquid water at operating point = water_condensed + free_water
        (stb/MMscf | sm3/sm3). This is the water that must be treated with
        inhibitor to prevent hydrate formation.

    Attributes — Inhibitor Injection Rate
    --------------------------------------
    Injection rate is based on the total liquid water at the operating point
    and the required inhibitor concentration (wt% in the aqueous phase).

    inhibitor_mass_rate : float
        Required inhibitor mass injection rate (lb/MMscf | kg/sm3 if metric).
        Zero when outside hydrate window, no inhibitor, or no liquid water.
    inhibitor_vol_rate : float
        Required inhibitor volume injection rate (gal/MMscf | L/sm3 if metric).
        Zero when outside hydrate window, no inhibitor, or no liquid water.
    """
    hft: float
    hfp: float
    subcooling: float
    in_hydrate_window: bool
    inhibited_hft: float
    inhibitor_depression: float
    required_inhibitor_wt_pct: float
    max_inhibitor_wt_pct: float
    inhibitor_underdosed: bool
    water_vaporized_res: float
    water_vaporized_op: float
    water_condensed: float
    free_water: float
    total_liquid_water: float
    inhibitor_mass_rate: float
    inhibitor_vol_rate: float


def _motiee_hft(p_psia, sg):
    """Motiee (1991) hydrate formation temperature.

    Published form uses degC and kPa with log10.
    T(degC) = -283.24469 + 78.99667*log10(P_kPa) - 5.352544*log10(P_kPa)^2
              + 349.473877*gamma - 150.854675*gamma^2 - 27.604065*gamma*log10(P_kPa)

    Reference: Motiee, M. (1991). Hydrocarbon Processing 70, pp 98-99.
    """
    p_kpa = p_psia * _PSI_TO_KPA
    log_p = math.log10(p_kpa)
    t_c = (_MOTIEE[0]
           + _MOTIEE[1] * log_p
           + _MOTIEE[2] * log_p * log_p
           + _MOTIEE[3] * sg
           + _MOTIEE[4] * sg * sg
           + _MOTIEE[5] * sg * log_p)
    return t_c * 9.0 / 5.0 + 32.0  # degC -> degF


def _towler_mokhatab_hft(p_psia, sg):
    """Towler & Mokhatab (2005) hydrate formation temperature.

    T(degF) = 13.47*ln(P_psia) + 34.27*ln(gamma) - 1.675*ln(P_psia)*ln(gamma) - 20.35

    Reference: Towler, B.F. & Mokhatab, S. (2005). Hydrocarbon Processing 84, pp 61-62.
    """
    ln_p = math.log(p_psia)
    ln_sg = math.log(sg)
    return _TOWLER[0] * ln_p + _TOWLER[1] * ln_sg + _TOWLER[2] * ln_p * ln_sg + _TOWLER[3]


def _hydrate_formation_press(degf_target, sg, hft_fn):
    """Invert HFT correlation to find hydrate formation pressure via bisection.

    Uses Motiee for inversion (consistent with ResToolbox3).
    Returns pressure in psia, or NaN if target T is outside correlation range.
    """
    p_lo = _HFP_P_LO
    p_hi = _HFP_P_HI

    t_lo = hft_fn(p_lo, sg)
    t_hi = hft_fn(p_hi, sg)

    # Check that solution exists within bounds
    if degf_target < t_lo or degf_target > t_hi:
        return float('nan')

    for _ in range(100):
        p_mid = (p_lo + p_hi) / 2.0
        t_mid = hft_fn(p_mid, sg)
        if abs(t_mid - degf_target) < 0.001:
            break
        if t_mid < degf_target:
            p_lo = p_mid
        else:
            p_hi = p_mid

    return (p_lo + p_hi) / 2.0


def _ostergaard_depression(wt_pct, inh):
    """Østergaard et al. (2005) temperature depression in degC.

    delta_T(degC) = C1*w + C2*w^2 + C3*w^3, where w = wt% (0-100 scale).

    Reference: Østergaard, K.K. et al. (2005). J. Pet. Sci. Eng. 48, pp 70-80.
    """
    c1, c2, c3 = _OSTERGAARD[inh]
    return c1 * wt_pct + c2 * wt_pct**2 + c3 * wt_pct**3


def _required_concentration(depression_degc, inh):
    """Newton-Raphson inversion of Østergaard cubic to find required wt%.

    Returns wt% (0-100 scale). Returns 0 if depression <= 0.
    """
    if depression_degc <= 0:
        return 0.0

    c1, c2, c3 = _OSTERGAARD[inh]

    # Initial guess from linear term
    w = depression_degc / c1 if c1 > 0 else 20.0

    for _ in range(50):
        f = c1 * w + c2 * w**2 + c3 * w**3 - depression_degc
        fp = c1 + 2 * c2 * w + 3 * c3 * w**2
        if abs(fp) < 1e-15:
            break
        w_new = w - f / fp
        if w_new < 0:
            w_new = w / 2.0
        if abs(w_new - w) < 1e-6:
            w = w_new
            break
        w = w_new

    return max(float(w), 0.0)


def gas_hydrate(
    p: float,
    degf: float,
    sg: float,
    hydmethod: str = 'TOWLER',
    inhibitor_type: str = None,
    inhibitor_wt_pct: float = 0,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    p_res: float = None,
    degf_res: float = None,
    additional_water: float = 0,
    metric: bool = False,
) -> HydrateResult:
    """ Returns gas hydrate formation prediction, water balance, and inhibitor calculations.

        Hydrate assessment (HFT, HFP, subcooling, inhibitor) is evaluated at
        the operating point (p, degf) — typically the wellhead or coldest point
        in the production system. Water balance is evaluated between reservoir
        conditions (p_res, degf_res) and the operating point to determine how
        much water condenses from vapor, how much was always liquid (free water),
        and the total liquid water that must be treated with inhibitor.

        p: Operating pressure at hydrate assessment point, e.g. wellhead
           (psia | barsa if metric=True)
        degf: Operating temperature at hydrate assessment point
              (degF | degC if metric=True)
        sg: Gas specific gravity (air = 1.0)
        hydmethod: Hydrate formation correlation.
                   'TOWLER': Towler & Mokhatab (2005)
                   'MOTIEE': Motiee (1991)
                   Defaults to 'TOWLER'
        inhibitor_type: Thermodynamic hydrate inhibitor type (optional).
                        'MEOH' (Methanol), 'MEG' (Monoethylene Glycol),
                        'DEG' (Diethylene Glycol), 'TEG' (Triethylene Glycol),
                        'ETOH' (Ethanol). None = no inhibitor
        inhibitor_wt_pct: Weight percent of inhibitor in aqueous phase (water
                          + inhibitor, 0-100). Defaults to 0
        co2: CO2 mole fraction (0-1). For composition-aware water content via
             SoreideWhitson. Defaults to 0
        h2s: H2S mole fraction (0-1). Defaults to 0
        n2: N2 mole fraction (0-1). Defaults to 0
        h2: H2 mole fraction (0-1). Defaults to 0
        p_res: Reservoir pressure where gas was last in equilibrium with water
               (psia | barsa if metric=True). Determines how much water the gas
               carries as vapor from the reservoir. If None, uses p (operating
               pressure — no condensation). Defaults to None
        degf_res: Reservoir temperature where gas was last in equilibrium with
                  water (degF | degC if metric=True). Determines how much water
                  the gas carries as vapor from the reservoir. If None, uses degf
                  (operating temperature — no condensation). Defaults to None
        additional_water: Free liquid water entrained in the gas stream from the
                          reservoir, e.g. from mobile formation water
                          (stb/MMscf | sm3/sm3 if metric). This water was never
                          vaporized — it travels with the gas as liquid. Added to
                          condensed water for inhibitor dosing. Defaults to 0
        metric: If True, input/output in Eclipse METRIC units (barsa, degC).
                Defaults to False (FIELD: psia, degF)

    Returns a HydrateResult dataclass. See HydrateResult docstring for full
    field descriptions.
    """
    # Convert metric inputs to oilfield
    if metric:
        p_psia = p * BAR_TO_PSI
        degf_of = degc_to_degf(degf)
        additional_water_stb = additional_water * SM3_PER_SM3_TO_STB_PER_MMSCF
        p_res_psia = p_res * BAR_TO_PSI if p_res is not None else None
        degf_res_of = degc_to_degf(degf_res) if degf_res is not None else None
    else:
        p_psia = p
        degf_of = degf
        additional_water_stb = additional_water
        p_res_psia = p_res
        degf_res_of = degf_res

    # Water content P,T: use reservoir conditions if provided, else operating
    p_wc = p_res_psia if p_res_psia is not None else p_psia
    degf_wc = degf_res_of if degf_res_of is not None else degf_of

    # Validate inputs
    validate_pe_inputs(p=p_psia, degf=degf_of, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    if p_wc != p_psia:
        validate_pe_inputs(p=p_wc)
    if inhibitor_wt_pct < 0 or inhibitor_wt_pct >= 100:
        raise ValueError("Inhibitor wt% must be in range [0, 100)")
    if additional_water < 0:
        raise ValueError("additional_water must be >= 0")

    # Resolve hydrate method
    hydmethod = validate_methods(["hydmethod"], [hydmethod])

    # Select HFT function
    if hydmethod == hyd_method.MOTIEE:
        hft_fn = _motiee_hft
    else:
        hft_fn = _towler_mokhatab_hft

    # Compute HFT and HFP
    hft_degf = hft_fn(p_psia, sg)
    hfp_psia = _hydrate_formation_press(degf_of, sg, hft_fn)

    # Subcooling and hydrate window
    subcooling_degf = hft_degf - degf_of
    in_window = degf_of < hft_degf

    # --- Water balance ---
    # Compute vaporized water at both reservoir and operating conditions.
    # Condensed water = what dropped out of vapor between the two points.
    # Free water = liquid water entrained from reservoir (user input).
    # Total liquid = condensed + free = what needs inhibitor treatment.
    has_composition = co2 > 0 or h2s > 0 or n2 > 0 or h2 > 0

    def _water_content_at(p_eval, degf_eval):
        """Compute equilibrium vaporized water content at given P,T (stb/MMscf)."""
        if has_composition:
            from pyrestoolbox.brine import SoreideWhitson
            sw = SoreideWhitson(
                pres=p_eval, temp=degf_eval, ppm=0,
                y_CO2=co2, y_H2S=h2s, y_N2=n2, y_H2=h2,
                sg=sg, metric=False,
            )
            return float(sw.water_content['stb_mmscf'])
        else:
            return float(gas_water_content(p=p_eval, degf=degf_eval, salinity=0, metric=False))

    # Vaporized water at operating point (always needed)
    wc_op = _water_content_at(p_psia, degf_of)

    # Vaporized water at reservoir (= operating if no reservoir P,T given)
    if p_res_psia is not None or degf_res_of is not None:
        wc_res = _water_content_at(p_wc, degf_wc)
    else:
        wc_res = wc_op  # no reservoir specified → no condensation

    # Condensed water: what dropped out of the gas between reservoir and operating
    condensed = max(wc_res - wc_op, 0.0)

    # Free water: liquid water entrained from reservoir (user input)
    free_water_stb = additional_water_stb

    # Total liquid water at operating point: this is what needs inhibitor
    total_liquid = condensed + free_water_stb

    # --- Inhibitor calculations ---
    inhibited_hft_degf = float('nan')
    depression_degf = 0.0
    required_wt_pct = 0.0
    max_wt_pct = 0.0
    underdosed = False
    inh_mass_rate = 0.0
    inh_vol_rate = 0.0

    if inhibitor_type is not None:
        # Resolve inhibitor enum
        inh = validate_methods(["inhibitor"], [inhibitor_type])
        max_wt_pct = _MAX_WT_PCT[inh]

        if inhibitor_wt_pct > 0:
            # Østergaard depression (in degC), convert to degF delta
            depression_degc = _ostergaard_depression(inhibitor_wt_pct, inh)
            depression_degf = depression_degc * 9.0 / 5.0
            inhibited_hft_degf = hft_degf - depression_degf
        else:
            inhibited_hft_degf = hft_degf
            depression_degf = 0.0

        # Required concentration to bring HFT below operating T (with capping)
        if degf_of < hft_degf:
            needed_depression_degc = (hft_degf - degf_of) * 5.0 / 9.0
            raw_wt_pct = _required_concentration(needed_depression_degc, inh)
            if raw_wt_pct > max_wt_pct:
                required_wt_pct = max_wt_pct
                underdosed = True
            else:
                required_wt_pct = raw_wt_pct
                underdosed = False
        else:
            required_wt_pct = 0.0

        # --- Injection rate ---
        # Inhibitor treats the total liquid water at operating conditions
        if in_window and required_wt_pct > 0 and total_liquid > 0:
            w_frac = required_wt_pct / 100.0
            liquid_mass_lb = total_liquid * _WATER_LB_PER_STB  # lb per MMscf gas
            inh_mass_rate = liquid_mass_lb * w_frac / (1.0 - w_frac)  # lb/MMscf
            density_lb_per_gal = _INHIBITOR_DENSITY[inh] * _GCM3_TO_LB_PER_GAL
            inh_vol_rate = inh_mass_rate / density_lb_per_gal  # gal/MMscf

    # Convert outputs if metric
    _wc_conv = STB_PER_MMSCF_TO_SM3_PER_SM3
    if metric:
        hft_out = degf_to_degc(hft_degf)
        hfp_out = hfp_psia * PSI_TO_BAR if not math.isnan(hfp_psia) else float('nan')
        subcooling_out = subcooling_degf * 5.0 / 9.0
        depression_out = depression_degf * 5.0 / 9.0
        inhibited_hft_out = degf_to_degc(inhibited_hft_degf) if not math.isnan(inhibited_hft_degf) else float('nan')
        wc_res_out = wc_res * _wc_conv
        wc_op_out = wc_op * _wc_conv
        condensed_out = condensed * _wc_conv
        free_water_out = free_water_stb * _wc_conv
        total_liquid_out = total_liquid * _wc_conv
        inh_mass_rate_out = inh_mass_rate * LB_PER_MMSCF_TO_KG_PER_SM3
        inh_vol_rate_out = inh_vol_rate * GAL_PER_MMSCF_TO_L_PER_SM3
    else:
        hft_out = hft_degf
        hfp_out = hfp_psia
        subcooling_out = subcooling_degf
        depression_out = depression_degf
        inhibited_hft_out = inhibited_hft_degf
        wc_res_out = wc_res
        wc_op_out = wc_op
        condensed_out = condensed
        free_water_out = free_water_stb
        total_liquid_out = total_liquid
        inh_mass_rate_out = inh_mass_rate
        inh_vol_rate_out = inh_vol_rate

    return HydrateResult(
        hft=hft_out,
        hfp=hfp_out,
        subcooling=subcooling_out,
        in_hydrate_window=in_window,
        inhibited_hft=inhibited_hft_out,
        inhibitor_depression=depression_out,
        required_inhibitor_wt_pct=required_wt_pct,
        max_inhibitor_wt_pct=max_wt_pct,
        inhibitor_underdosed=underdosed,
        water_vaporized_res=wc_res_out,
        water_vaporized_op=wc_op_out,
        water_condensed=condensed_out,
        free_water=free_water_out,
        total_liquid_water=total_liquid_out,
        inhibitor_mass_rate=inh_mass_rate_out,
        inhibitor_vol_rate=inh_vol_rate_out,
    )
