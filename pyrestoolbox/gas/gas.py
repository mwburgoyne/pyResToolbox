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

Classes
-------
GasPVT          Convenience wrapper storing gas composition & method choices
"""

__all__ = [
    'gas_z', 'gas_ug', 'gas_bg', 'gas_cg', 'gas_den', 'gas_sg',
    'gas_tc_pc', 'gas_dmp', 'gas_ponz2p', 'gas_grad2sg', 'gas_fws_sg',
    'gas_water_content', 'gas_rate_radial', 'gas_rate_linear', 'darcy_gas',
    'GasPVT',
    # Enum classes re-exported for convenience (gas.z_method.DAK, etc.)
    'z_method', 'c_method',
]

import numpy as np
import numpy.typing as npt
from typing import Tuple

import pandas as pd
from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic
from pyrestoolbox.shared_fns import convert_to_numpy, process_output, check_2_inputs, bisect_solve, validate_pe_inputs
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.constants import (R, psc, tsc, degF2R, tscr, scf_per_mol, CUFTperBBL, WDEN, MW_CO2, MW_H2S, MW_N2, MW_AIR, MW_H2,
    BAR_TO_PSI, PSI_TO_BAR, degc_to_degf, degf_to_degc,
    M_TO_FT, FT_TO_M, MM_TO_IN, IN_TO_MM, SQM_TO_SQFT, SQFT_TO_SQM,
    LBCUFT_TO_KGM3, KGM3_TO_LBCUFT, INVPSI_TO_INVBAR, INVBAR_TO_INVPSI,
    PSI2CP_TO_BAR2CP, BAR2CP_TO_PSI2CP, BARM_TO_PSIFT, PSIFT_TO_BARM,
    MSCF_TO_SM3, SM3_TO_MSCF, STB_PER_MMSCF_TO_SM3_PER_SM3,
    SM3_PER_SM3_TO_STB_PER_MSCF, D_PER_SM3_TO_D_PER_MSCF, D_PER_MSCF_TO_D_PER_SM3)

# Precomputed Gauss-Legendre nodes/weights for pseudopressure integration
_GL7_NODES, _GL7_WEIGHTS = np.polynomial.legendre.leggauss(7)
_GL10_NODES, _GL10_WEIGHTS = np.polynomial.legendre.leggauss(10)

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
    if h2 > 0:
        cmethod = 'BNS'
        zmethod = 'BNS'
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])
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
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, m, sm3/d). Defaults to False (FIELD)
    """
    if metric:
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        h = np.asarray(h) * M_TO_FT if not isinstance(h, (int, float)) else h * M_TO_FT
        degf = degc_to_degf(degf)
        r_w = r_w * M_TO_FT
        r_ext = r_ext * M_TO_FT
        if D > 0:
            D = D * D_PER_SM3_TO_D_PER_MSCF  # day/sm3 -> day/Mscf
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    if r_w <= 0:
        raise ValueError("Wellbore radius r_w must be positive")
    if r_ext <= r_w:
        raise ValueError("External radius r_ext must be greater than wellbore radius r_w")

    k, h, pr, pwf = np.asarray(k), np.asarray(h), np.asarray(pr), np.asarray(pwf)
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
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, m, m2, sm3/d). Defaults to False (FIELD)
    """
    if metric:
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        area = np.asarray(area) * SQM_TO_SQFT if not isinstance(area, (int, float)) else area * SQM_TO_SQFT
        degf = degc_to_degf(degf)
        length = length * M_TO_FT
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    if length <= 0:
        raise ValueError("Flow length must be positive")

    k, area, pr, pwf = np.asarray(k), np.asarray(area), np.asarray(pr), np.asarray(pwf)
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
    tr = degf + degF2R
    if radial:
        a = k * h * delta_mp
        b = 1422 * tr
        c = np.log(l2 / l1) - 0.75 + S
        if D > 1e-9:  # Solve analytically for rate with non-Darcy factor by rearranging into root of a quadratic equation.
            return (np.sqrt(4 * a * b * D + (b * b * c * c)) - (b * c)) / (2 * b * D)
    else:
        a = k * h * l1 * delta_mp
        b = 1422 * tr
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
    """ Returns a Tuple of critical temperature (deg R) and critical pressure (psia) for hydrocarbon gas
        For SUT and PMC, this returns an equivalent set of critical parameters for the mixture
        For BUR, this returns critical parameters for the pure hydrocarbon gas alone
        sg: Specific gravity of reservoir gas (relative to air)
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        cmethod: 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, 
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BNS' for Burgoyne, Nielsen and Stanko method (2025). If h2 > 0, then 'BNS' will be used
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified
        metric: If True, input/output in Eclipse METRIC units (K, barsa). Defaults to False (FIELD)
    """
    if metric:
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    if tc * pc > 0:  # Critical properties have both been user specified
        if metric:
            return (tc / 1.8, pc * PSI_TO_BAR)  # deg R -> K, psia -> barsa
        return (tc, pc)
    
    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
    
    cmethod = validate_methods(["cmethod"], [cmethod])

    if cmethod.name == "PMC":  # Piper, McCain & Corredor (1999)
        y = np.array([0, h2s, co2, n2])
        alpha = np.array(
            [0.11582, -0.4582, -0.90348, -0.66026, 0.70729, -0.099397]
        )
        beta = np.array(
            [3.8216, -0.06534, -0.42113, -0.91249, 17.438, -3.2191]
        )
        tci = np.array([0, 672.35, 547.58, 239.26])
        pci = np.array([0, 1306.0, 1071.0, 507.5])
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
        if hc_frac <= 0:
            raise ValueError("SUT method requires hydrocarbon fraction > 0 (n2 + co2 + h2s must be < 1.0)")
        sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / MW_AIR) / hc_frac  # Eq 3.53
        ppc_hc = 756.8 - 131.0 * sg_hc - 3.6 * sg_hc ** 2  # Eq 3.47b
        tpc_hc = 169.2 + 349.5 * sg_hc - 74.0 * sg_hc ** 2  # Eq 3.47a

        # Wichert & Aziz non-hydrocarbon corrections from monograph
        eps = 120 * ((co2 + h2s) ** 0.9 - (co2 + h2s) ** 1.6) + 15 * (
            h2s ** 0.5 - h2s ** 4
        )  # Eq 3.52c
        ppc_star = (
            (1 - n2 - co2 - h2s) * ppc_hc
            + n2 * 507.5
            + co2 * 1071.0
            + h2s * 1306.0
        )  # Eq 3.54a
        tpc_star = (
            (1 - n2 - co2 - h2s) * tpc_hc
            + n2 * 239.26
            + co2 * 547.58
            + h2s * 672.35
        )  # Eq 3.54b
        tpc = tpc_star - eps  # Eq 3.52a
        ppc = (
            ppc_star * (tpc_star - eps) / (tpc_star + h2s * (1 - h2s) * eps)
        )  # Eq. 3.52b


    elif cmethod.name in ("BUR", "BNS"):
        def tc_ag(x):
            a, b, c = 2695.14765, 274.341701, 343.008
            return a * x / (b + x) + c
        
        def tc_gc(x):
            a, b, c = 1098.10948, 101.529237, 343.008
            return a * x / (b + x) + c
        
        def pc_fn(x, vc_slope, tc_p):
            vc_on_zc = vc_slope * x + 5.518525872412144
            return R * tc_p / vc_on_zc
        
        def pseudo_critical(sg_hc, AG = False):
            """
            Calculates pseudo-critical temperature and pressure for the pseudo-hydrocarbon component.
            - Uses custom linear fits for gas condensates (AG = False) and associated gas (AG = True) (Burgoyne, 2025).
            - These relations are derived from property fitting and may differ from literature.
            """
            x = max(0, MW_AIR * sg_hc - 16.0425)
            if AG:
                tpc_hc = tc_ag(x)
                vc_slope = 0.177497835
            else:
                tpc_hc = tc_gc(x)
                vc_slope = 0.170931432
            ppc_hc = pc_fn(x, vc_slope, tpc_hc)
            return tpc_hc, ppc_hc   
            
        if co2 + h2s + n2 + h2 < 1.0: # If not 100% Inerts, then calculate hydrocarbon MW
            hydrocarbon_specific_gravity = (sg - (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2) / MW_AIR) / (1 - co2 - h2s - n2 - h2)
        else:
            hydrocarbon_specific_gravity = 0.75 # Use default value if 100% inerts to avoid numerical problems
        hydrocarbon_specific_gravity = np.max([0.553779772, hydrocarbon_specific_gravity])  # Methane is lower limit
        
        hydrocarbon_molecular_weight = hydrocarbon_specific_gravity * MW_AIR
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

def _halley_cubic_vec(c2, c1, c0, max_iter=50, tol=1e-12):
    """Vectorized Halley solver: solve Z^3+c2*Z^2+c1*Z+c0=0 for max root (vapor Z).
    c2, c1, c0 are 1D arrays of length N. Returns 1D array of Z values.
    Falls back to _cardano_cubic for any non-converged elements."""
    N = len(c2)
    Z = -c2 / 3.0
    f0 = Z**3 + c2 * Z**2 + c1 * Z + c0
    Z = np.where(f0 < 0, Z + 1.0, Z)

    for _ in range(max_iter):
        f = Z**3 + c2 * Z**2 + c1 * Z + c0
        fp = 3.0 * Z**2 + 2.0 * c2 * Z + c1
        fpp = 6.0 * Z + 2.0 * c2
        # Protect against zero derivatives
        safe_fp = np.where(np.abs(fp) < 1e-30, 1e-30, fp)
        dZ = f / safe_fp
        denom = safe_fp - 0.5 * dZ * fpp
        denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
        dZ = f / denom
        Z -= dZ
        if np.max(np.abs(dZ)) < tol:
            break

    # Check residuals and fall back to Cardano for any bad elements
    f = Z**3 + c2 * Z**2 + c1 * Z + c0
    bad = np.abs(f) > 1e-6
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
        tc: Critical gas temperature (deg R | K). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia | barsa). Uses cmethod correlation if not specified
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, K). Defaults to False (FIELD)
    """
    if metric:
        p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    validate_pe_inputs(p=p, degf=degf, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)

    tolerance = 1e-6
    p, is_list = convert_to_numpy(p)

    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS' 
    elif cmethod == 'BNS':
        zmethod = 'BNS'
    elif zmethod == 'BNS':
        cmethod='BNS'    
        
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    
    if n2 + co2 + h2s + h2 < 1:
        sg_hc = (sg - (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2) / MW_AIR) / (1 - co2 - h2s - n2 - h2)
    else:
        sg_hc = 0.75 # Irrelevant, since hydrocarbon fraction = 0
    sg_hc = max(sg_hc, 0.553779772) # Methane is lower limit
    mw_hc = sg_hc * MW_AIR
        
    tc, pc = gas_tc_pc(sg, co2, h2s, n2, h2, cmethod.name, tc, pc)
    tr = (degf + degF2R) / tc
    pprs = np.array(p/pc)
            
    def zdak(pprs, tr):
        # DAK from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        # Vectorized Newton-Raphson on reduced density rhor
        A1, A2, A3, A4, A5 = 0.3265, -1.0700, -0.5339, 0.01569, -0.05165
        A6, A7, A8, A9, A10, A11 = 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210

        R1 = A1 + A2 / tr + A3 / tr**3 + A4 / tr**4 + A5 / tr**5
        R2 = 0.27 * pprs / tr  # Array
        R3 = A6 + A7 / tr + A8 / tr**2
        R4 = A9 * (A7 / tr + A8 / tr**2)
        R5 = A10 / tr**3

        rhor = 0.27 * pprs / tr  # Initial guess (array)
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

        zout = 0.27 * pprs / (rhor * tr)
        return process_output(zout, is_list)

    # Hall & Yarborough — Vectorized Newton-Raphson
    def z_hy(pprs, tr):

        tpr_inv = 1/tr  # Reciprocal reduced temperature
        t2 = tpr_inv ** 2
        a = 0.06125 * tpr_inv * np.exp(-1.2 * (1 - tpr_inv) ** 2)
        b = tpr_inv * (14.76 - 9.76 * tpr_inv + 4.58 * t2)
        c = tpr_inv * (90.7 - 242.2 * tpr_inv + 42.4 * t2)
        D = 2.18 + 2.82 * tpr_inv

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
        a = [0, 256.41675, 7.18202, -178.5725, 182.98704, -40.74427, 2.24427, 47.44825, 5.2852, -0.14914, 271.50446, 16.2694, -121.51728, 167.71477, -81.73093, 20.36191, -2.1177, 124.64444, -6.74331, 0.20897, -0.00314]
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
        Z_raw = _halley_cubic_vec(c2, c1_coeff, c0)    # (N,)

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
        
    zfuncs = {"DAK": zdak, "HY": z_hy, "WYW": z_wyw, "BUR": z_bur}

    if zmethod.name in ('BNS', 'BUR'):
        return zfuncs["BUR"](p, degf)
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
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
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

    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS'
    elif cmethod == 'BNS':
        zmethod = 'BNS'
    elif zmethod == 'BNS':
        cmethod='BNS'

    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    t = degf + degF2R
    m = MW_AIR * sg

    if not zee_provided or not check_2_inputs(zee, p): # Need to calculate Z-Factors if not same length / type as p
         zee = gas_z(p, sg, degf, zmethod, cmethod, co2, h2s, n2, h2, tc, pc)
    
    rho = m * p / (t * zee * R * 62.37)
    
    if zmethod.name not in ('BNS', 'BUR'):
        b = 3.448 + (986.4 / t) + (0.01009 * m)  # 2.16
        c = 2.447 - (0.2224 * b)  # 2.17
        a = ((9.379 + (0.01607 * m)) * np.power(t, 1.5) / (209.2 + (19.26 * m) + t))  # 2.15
        ug = process_output(a * 0.0001 * np.exp(b * np.power(rho, c)), is_list)  # 2.14
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

        sg_hc = max(sg_hc, 0.553779772)
        hc_gas_mw = sg_hc * MW_AIR

        mws_lbc[-1] = hc_gas_mw
        tcs_lbc[-1], pcs_lbc[-1] = gas_tc_pc(hc_gas_mw / MW_AIR, cmethod='BNS')
        VCVIS_lbc[-1] = 0.0576710 * (hc_gas_mw - 16.0425) + 1.44383

        # Vectorized Stiel-Thodos
        Tr = degR / tcs_lbc
        Tc_K = tcs_lbc * 5.0 / 9.0
        Pc_atm = pcs_lbc / 14.696
        eta_st = Tc_K**(1.0/6.0) / (mws_lbc**0.5 * Pc_atm**(2.0/3.0))
        ui_low = 34e-5 * Tr**0.94 / eta_st
        ui_high = 17.78e-5 * np.maximum(4.58 * Tr - 1.67, 1e-30)**(5.0/8.0) / eta_st
        ui = np.where(Tr <= 1.5, ui_low, ui_high)

        # Herning-Zippener dilute gas mixture viscosity
        sqrt_mws = np.sqrt(mws_lbc)
        u0_val = np.sum(zi * ui * sqrt_mws) / np.sum(zi * sqrt_mws)

        # LBC mixture parameters
        a_lbc = np.array([0.1023, 0.023364, 0.058533, -0.0392852, 0.00926279])
        rhoc = 1.0 / np.sum(VCVIS_lbc * zi)
        eta_mix = np.abs(np.sum(zi * Tc_K))**(1.0/6.0) / (np.abs(np.sum(zi * mws_lbc))**0.5 * np.abs(np.sum(zi * Pc_atm))**(2.0/3.0))

        # Vectorized over pressures
        zee_arr, _ = convert_to_numpy(zee)
        rhor = p / (zee_arr * R * degR * rhoc)
        lhs = a_lbc[0] + rhor * (a_lbc[1] + rhor * (a_lbc[2] + rhor * (a_lbc[3] + rhor * a_lbc[4])))
        ug = process_output((lhs**4 - 1e-4) / eta_mix + u0_val, is_list)
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
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
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
    if metric:
        p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS'

    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    p, is_list = convert_to_numpy(p)
    tc, pc = gas_tc_pc(sg=sg, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc, cmethod=cmethod)

    pr = p / pc
    degR = (degf + degF2R)
    tr = degR / tc
    zee1 = gas_z(p=p, sg=sg, degf=degf, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2=h2, tc=tc, pc=pc)
    zee2 = gas_z(p=p+1, sg=sg, degf=degf,  zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2=h2, tc=tc, pc=pc)

    vol1 = zee1*R*degR/p
    vol2 = zee2*R*degR/(p+1)

    result = process_output((vol1 - vol2)/((vol1 + vol2)/2), is_list) # 1/psi
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
          tc: Critical gas temperature (deg R | K). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia | barsa). Calculates using cmethod if not specified
          metric: If True, input/output in Eclipse METRIC units (barsa, degC). Defaults to False (FIELD)
    """
    if metric:
        p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    p, is_list = convert_to_numpy(p)
    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS'

    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

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
          tc: Critical gas temperature (deg R | K). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia | barsa). Calculates using cmethod if not specified
          metric: If True, input/output in Eclipse METRIC units (barsa, degC, kg/m3). Defaults to False (FIELD)
    """
    if metric:
        p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if tc > 0:
            tc = tc * 1.8  # K to deg R
        if pc > 0:
            pc = pc * BAR_TO_PSI
    p, is_list = convert_to_numpy(p)

    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS'

    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

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
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
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

    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS'

    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    def PonZ2P_err(args, p):
        ponz, sg, degf, zmethod, cmethod, tc, pc, co2, h2s, n2, h2 = args
        zee = zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)
        return (p - (ponz * zee)) / p

    poverz, is_list = convert_to_numpy(poverz)

    p = []
    for ponz in poverz:
        args = (ponz, sg, degf, zmethod, cmethod, tc, pc, co2, h2s, n2, h2)
        p.append(bisect_solve(args, PonZ2P_err, ponz * 0.1, ponz * 5, rtol))

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
        Calculated through iterative solution method. Will fail if gas SG is below 0.55, or greater than 3.0.

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
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
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

    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS' 
        
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    args = (grad, p, zmethod, cmethod, tc, pc, co2, h2s, n2, h2)
    return bisect_solve(args, grad_err, 0.06958923023817742, 3.0, rtol)

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
        tc: Critical gas temperature (deg R | K). Calculates using cmethod if not specified
        pc: Critical gas pressure (psia | barsa). Calculates using cmethod if not specified
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
    if p1 == p2:
        return 0

    if h2 > 0:
        cmethod = 'BNS' # The BNS PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BNS'

    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

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
    cond_mw = 240 - 2.22 * api_st  # lb/lb-moles (Standing correlation)
    cond_moles = cond_mass / cond_mw  # lb-moles
    fws_gas_mass = cond_mass + surface_gas_mass  # lb
    fws_gas_moles = cond_moles + surface_gas_moles  # lb-moles
    fws_gas_mw = fws_gas_mass / fws_gas_moles  # lb/lb-moles
    return fws_gas_mw / MW_AIR

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
        metric: If True, methods accept/return Eclipse METRIC units (barsa, degC, kg/m3). Defaults to False
    """
    def __init__(self, sg=0.75, co2=0, h2s=0, n2=0, h2=0,
                 zmethod='DAK', cmethod='PMC', metric=False):
        self.sg = sg
        self.co2 = co2
        self.h2s = h2s
        self.n2 = n2
        self.h2 = h2
        self.metric = metric
        if h2 > 0:
            zmethod = 'BNS'
            cmethod = 'BNS'
        self.zmethod, self.cmethod = validate_methods(
            ["zmethod", "cmethod"], [zmethod, cmethod])
        # tc/pc always stored in oilfield units internally
        self.tc, self.pc = gas_tc_pc(sg, co2, h2s, n2, h2, self.cmethod.name)

    def _convert_inputs(self, p, degf):
        """Convert metric inputs to oilfield for internal calculations."""
        if self.metric:
            p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI
            degf = degc_to_degf(degf)
        return p, degf

    def z(self, p, degf):
        """ Returns Z-factor at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        return gas_z(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                     cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                     n2=self.n2, h2=self.h2, tc=self.tc, pc=self.pc)

    def viscosity(self, p, degf):
        """ Returns gas viscosity (cP) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        return gas_ug(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                      cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                      n2=self.n2, h2=self.h2, tc=self.tc, pc=self.pc)

    def density(self, p, degf):
        """ Returns gas density (lb/cuft | kg/m3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        result = gas_den(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                         cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                         n2=self.n2, h2=self.h2, tc=self.tc, pc=self.pc)
        if self.metric:
            return result * LBCUFT_TO_KGM3
        return result

    def bg(self, p, degf):
        """ Returns gas FVF Bg (rcf/scf | rm3/sm3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        p, degf = self._convert_inputs(p, degf)
        return gas_bg(p=p, sg=self.sg, degf=degf, zmethod=self.zmethod,
                      cmethod=self.cmethod, co2=self.co2, h2s=self.h2s,
                      n2=self.n2, h2=self.h2, tc=self.tc, pc=self.pc)


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
            47484
            * (
                np.exp(
                    69.103501
                    + (-13064.76 / (t + degF2R))
                    + (-7.3037 * np.log(t + degF2R))
                    + (0.0000012856 * ((t + degF2R) * (t + degF2R)))
                )
            )
            / (p)
            + (np.power(10, ((-3083.87 / (t + degF2R)) + 6.69449)))
        )
        * (1 - (0.00492 * salinity) - (0.00017672 * (salinity * salinity)))
        / 8.32
        / 42
    )
    if metric:
        return content * STB_PER_MMSCF_TO_SM3_PER_SM3  # stb/MMscf -> sm3/sm3
    return content
