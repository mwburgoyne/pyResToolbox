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

Oil PVT, flow rate, and black oil table calculations.

Functions
---------
oil_pbub            Bubble point pressure (Standing, Valko-McCain, Velarde)
oil_rs_bub          Solution GOR at bubble point
oil_rs              Solution GOR at any pressure
oil_bo              Oil formation volume factor (McCain, Standing)
oil_deno            Live oil density
oil_viso            Oil viscosity (saturated and undersaturated)
oil_co              Oil compressibility
oil_sg              Oil specific gravity from API
oil_api             API gravity from specific gravity
oil_ja_sg           Jacoby aromaticity SG
oil_twu_props       Twu critical property correlations
oil_rs_st           Standing Rs correlation
oil_rate_radial     Radial oil flow rate (STB/d)
oil_rate_linear     Linear oil flow rate (STB/d)
oil_harmonize_pb_rsb  Harmonize consistent Pb and Rsb
sg_evolved_gas      Evolved gas specific gravity
sg_st_gas           Stock-tank gas specific gravity
sgg_wt_avg          Weighted average gas SG from separator stages
check_sgs           Validate separator/stock-tank gas SG consistency
make_bot_og         Black oil table generation (backward-compatible wrapper)

Classes
-------
OilPVT              Convenience wrapper storing oil characterization & method choices
"""

__all__ = [
    'oil_pbub', 'oil_rs_bub', 'oil_rs', 'oil_bo', 'oil_deno', 'oil_viso',
    'oil_co', 'oil_sg', 'oil_api', 'oil_ja_sg', 'oil_twu_props', 'oil_rs_st',
    'oil_rate_radial', 'oil_rate_linear', 'oil_harmonize_pb_rsb',
    'sg_evolved_gas', 'sg_st_gas', 'sgg_wt_avg', 'check_sgs',
    'make_bot_og', 'OilPVT',
    # Enum classes re-exported for convenience (oil.pb_method.STAN, etc.)
    'pb_method', 'rs_method', 'bo_method', 'co_method',
]

import numpy as np
import numpy.typing as npt
import pandas as pd
from tabulate import tabulate
from typing import Tuple

from pyrestoolbox.constants import R, psc, tsc, degF2R, tscr, MW_AIR, scf_per_mol, CUFTperBBL, WDEN
from pyrestoolbox.constants import (BAR_TO_PSI, PSI_TO_BAR, degc_to_degf, degf_to_degc,
    M_TO_FT, FT_TO_M, SQM_TO_SQFT, LBCUFT_TO_KGM3, KGM3_TO_LBCUFT,
    INVPSI_TO_INVBAR, INVBAR_TO_INVPSI,
    SCF_PER_STB_TO_SM3_PER_SM3, SM3_PER_SM3_TO_SCF_PER_STB,
    STB_TO_SM3, SM3_TO_STB)
from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic
from pyrestoolbox.validate import validate_methods
import pyrestoolbox.gas as gas
import pyrestoolbox.brine as brine

def get_real_part(value):
    if isinstance(value, complex):
        return value.real
    else:
        return value

def oil_sg(api_value: float) -> float:
    """ Returns oil specific gravity given API value of oil
        api_value: API value
    """
    return 141.5 / (api_value + 131.5)

def oil_api(sg_value: float) -> float:
    """ Returns oil API given specific gravity value of oil
        sg_value: Specific gravity (relative to water)
    """
    return 141.5 / sg_value - 131.5
    
def oil_rate_radial(
    k: npt.ArrayLike,
    h: npt.ArrayLike,
    pr: npt.ArrayLike,
    pwf: npt.ArrayLike,
    r_w: float,
    r_ext: float,
    uo: float,
    bo: float,
    S: float = 0,
    vogel: bool = False,
    pb: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns liquid rate for radial flow (stb/day | sm3/day) using Darcy pseudo steady state equation
        k: Effective Permeability to flow (mD)
        h: Net flow height (ft | m)
        Pr: Reservoir pressure (psia | barsa)
        pwf: BHFP (psia | barsa)
        r_w: Wellbore Radius (ft | m)
        r_ext: External Reservoir Radius (ft | m)
        uo: Liquid viscosity (cP)
        bo: Liquid Formation Volume Factor (rb/stb | rm3/sm3)
        S: Wellbore Skin (Dimensionless). Defaults to zero if not specified
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
        pb: Bubble point pressure (psia | barsa). Defaults to zero if not specified. Not used unless Vogel option is invoked
        metric: If True, input/output in Eclipse METRIC units (barsa, m, sm3/d). Defaults to False (FIELD)
    """
    if metric:
        h = np.asarray(h) * M_TO_FT if not isinstance(h, (int, float)) else h * M_TO_FT
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        r_w = r_w * M_TO_FT
        r_ext = r_ext * M_TO_FT
        if pb > 0:
            pb = pb * BAR_TO_PSI

    k, h, pr, pwf = (
        np.asarray(k),
        np.asarray(h),
        np.asarray(pr),
        np.asarray(pwf),
    )

    if pwf.size > 1:
        pb = np.array([max(pb, pwf[i]) for i in range(pwf.size)])

    if r_w <= 0:
        raise ValueError("Wellbore radius r_w must be positive")
    if r_ext <= r_w:
        raise ValueError("External radius r_ext must be greater than wellbore radius r_w")
    if uo <= 0:
        raise ValueError("Viscosity uo must be positive")
    if bo <= 0:
        raise ValueError("Formation volume factor bo must be positive")

    J = (
        0.00708 * k * h / (uo * bo * (np.log(r_ext / r_w) + S - 0.75))
    )  # Productivity index
    if not vogel:
        qoil = J * (pr - pwf)
    else:
        if np.any(pb <= 0):
            raise ValueError("Bubble point pressure pb must be positive when vogel=True")
        qsat_max = J * pb / 1.8
        qusat = J * (pr - pb)
        qoil = (
            qsat_max * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2) + qusat
        )
    if metric:
        return qoil * STB_TO_SM3  # stb/d -> sm3/d
    return qoil

def oil_rate_linear(
    k: npt.ArrayLike,
    pr: npt.ArrayLike,
    pwf: npt.ArrayLike,
    area: npt.ArrayLike,
    length: float,
    uo: float,
    bo: float,
    vogel: bool = False,
    pb: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns liquid rate for linear flow (stb/day | sm3/day) using Darcy steady state equation
        k: Permeability (mD)
        Pr: Reservoir pressure (psia | barsa)
        pwf: BHFP (psia | barsa)
        area: Net cross sectional area perpendicular to direction of flow (ft2 | m2).
        length: Length over which flow takes place (ft | m)
        uo: Liquid viscosity (cP)
        bo: Liquid Formation Volume Factor (rb/stb | rm3/sm3)
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
        pb: Bubble point pressure (psia | barsa). Defaults to zero if not specified. Not used unless Vogel option is invoked
        metric: If True, input/output in Eclipse METRIC units (barsa, m, m2, sm3/d). Defaults to False (FIELD)
    """
    if metric:
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        area = np.asarray(area) * SQM_TO_SQFT if not isinstance(area, (int, float)) else area * SQM_TO_SQFT
        length = length * M_TO_FT
        if pb > 0:
            pb = pb * BAR_TO_PSI

    k, area, pr, pwf = (
        np.asarray(k),
        np.asarray(area),
        np.asarray(pr),
        np.asarray(pwf),
    )

    if length <= 0:
        raise ValueError("Flow length must be positive")
    if uo <= 0:
        raise ValueError("Viscosity uo must be positive")
    if bo <= 0:
        raise ValueError("Formation volume factor bo must be positive")

    J = (
        0.00708 * k * area / (2 * np.pi * uo * bo * length)
    )  # Productivity index
    if not vogel:
        qoil = J * (pr - pwf)
    else:
        if np.any(pb <= 0):
            raise ValueError("Bubble point pressure pb must be positive when vogel=True")
        qsat_max = J * pb / 1.8
        qusat = J * (pr - pb)
        qoil = (
            qsat_max * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2) + qusat
        )
    if metric:
        return qoil * STB_TO_SM3  # stb/d -> sm3/d
    return qoil

def oil_ja_sg(mw: float, ja: float) -> float:
    """ Returns liquid hydrocarbon specific gravity using Jacoby Aromaticity Factor relationship
        mw: Molecular weight of the liquid (g/gmole or lb/lb-mol)
        Ja: Varies between 0 (Paraffins) - 1 (Aromatic)
    """
    ja = min(1, ja)
    ja = max(0, ja)
    return 0.8468 - 15.8 / mw + ja * (0.2456 - 1.77 / mw)


def oil_twu_props(
    mw: float, ja: float = 0, sg: float = 0, damp: float = 1, metric: bool = False
) -> Tuple:
    """ Returns Tuple of sg, tb (degR | K), tc (DegR | K), pc (psia | barsa), Vc (ft3/lb-mol) using method from Twu (1984) correlations for petroleum liquids
        Modified with damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu)
        Returns sg, tb (R | K), tc (R | K), pc (psia | barsa), vc (ft3/lbmol)

        mw: Molecular weight of the liquid hydrocarbon (g/g.mol / lb/lb.mol)
        ja: jacoby Aromaticity Factor relationship. Varies between 0 (Paraffins) - 1 (Aromatic). Defaults to zero if undefined
        sg: Specific gravity of the liquid (fraction relative to water density). Will use jacoby method to estimate sg from mw if undefined.
        damp: damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu). Defaults to 1
        metric: If True, output Tc/Tb in Kelvin and Pc in barsa. Defaults to False (FIELD: Rankine, psia)
        Unless otherwise mentioned, all Twu equation references are from Whitson Monograph
    """
    if sg == 0:
        sg = oil_ja_sg(mw, ja)  # Use jacoby relationship to estimate sg if not specified
        #print('sg', sg)
        
    # Estimate boiling point 
    # Return boiling point (deg R) and Paraffin properties
    def Twu_tb(mw, sg, damp=1):
        Mp_guess = mw  # Guess for paraffinic mw
        tb, tcp, pcp, vcp, sgp = paraffin_props(Mp_guess)
        d_err = mw - M(tb, sgp, sg, Mp_guess, damp)
        n_iter = 0
        while abs(d_err / mw) > 0.0001:
            n_iter += 1
            Mp_guess += d_err
            tb, tcp, pcp, vcp, sgp = paraffin_props(Mp_guess)
            d_err = mw - M(tb, sgp, sg, Mp_guess, damp)
            if n_iter > 100:
                raise RuntimeError(f"Twu algorithm did not converge for mw={mw}, ja={ja}, sg={sg}, damp={damp}")
        return tb, Mp_guess, tcp, pcp, vcp, sgp

    # Return mw from modified Eq 5.78 to take into account damping
    def M(tb, sgp, sg, Mp, damp):
        absx = abs(0.012342 - 0.328086 / tb ** 0.5)  # Just above Eq 5.78
        dsgM = (
            np.exp(5 * (sgp - sg)) - 1
        )  # Modified Eq 5.78 to take into account damping
        fm = dsgM * (
            absx + (-0.0175691 + 0.193168 / tb ** 0.5) * dsgM
        )  # Just above Eq 5.78
        M = np.exp(
            np.log(Mp) * (1 + 8 * damp * fm / (1 - 2 * fm) ** 2)
        )  # Modified Eq 5.78 to take into account damping
        return M

    def Twu_tc(tb, sgp, sg):
        tcp = (
            tb
            * (
                0.533272
                + 0.000191017 * tb
                + 0.0000000779681 * tb ** 2
                - 2.84376e-11 * tb ** 3
                + 95.9468 / (0.01 * tb) ** 13
            )
            ** -1
        )  # Eq 5.67
        dsgT = np.exp(5 * (sgp - sg)) - 1  # Eq 5.75
        ft = dsgT * (
            (-0.362456 / tb ** 0.5)
            + (0.0398285 - (0.948125 / tb ** 0.5)) * dsgT
        )  # Eq 5.75
        tc = tcp * ((1 + 2 * ft) / (1 - 2 * ft)) ** 2  # Eq 5.75
        return tc

    def Twu_vc(tb, tcp, sg, sgp):
        alpha = 1 - tb / tcp  # Eq 5.72
        vcp = (
            1
            - (
                0.419869
                - 0.505839 * alpha
                - 1.56436 * alpha ** 3
                - 9481.7 * alpha ** 14
            )
        ) ** -8  # Eq 5.69
        dsgV = np.exp(4 * (sgp ** 2 - sg ** 2)) - 1  # Eq 5.76
        f_v = dsgV * (
            (0.46659 / tb ** 0.5) + (-0.182421 + (3.01721 / tb ** 0.5)) * dsgV
        )  # Eq 5.76
        vc = vcp * ((1 + 2 * f_v) / (1 - 2 * f_v)) ** 2  # Eq 5.76
        return vc

    def Twu_pc(tb, sgp, sg, pcp, tc, tcp, vc, vcp):
        dsgp = np.exp(0.5 * (sgp - sg)) - 1  # Eq 5.77
        fp = dsgp * (
            (2.53262 - 46.1955 / tb ** 0.5 - 0.00127885 * tb)
            + (-11.4277 + 252.14 / tb ** 0.5 + 0.00230533 * tb) * dsgp
        )  # Eq 5.77
        pc = (
            pcp * (tc / tcp) * (vcp / vc) * ((1 + 2 * fp) / (1 - 2 * fp)) ** 2
        )  # Eq 5.77
        return pc

    def paraffin_props(Mp):
        theta = np.log(Mp)  # Eq 5.73
        tb = (
            np.exp(
                5.71419
                + 2.71579 * theta
                - 0.28659 * theta ** 2
                - 39.8544 / theta
                - 0.122488 / theta ** 2
            )
            - 24.7522 * theta
            + 35.3155 * theta ** 2
        )  # Eq 5.71
        tcp = (
            tb
            * (
                0.533272
                + 0.000191017 * tb
                + 0.0000000779681 * tb ** 2
                - 2.84376e-11 * tb ** 3
                + 95.9468 / (0.01 * tb) ** 13
            )
            ** -1
        )  # Eq. 5.67
        alpha = 1 - tb / tcp  # Eq 5.72
        pcp = (
            3.83354
            + 1.19629 * alpha ** 0.5
            + 34.8888 * alpha
            + 36.1952 * alpha ** 2
            + 104.193 * alpha ** 4
        ) ** 2  # Eq 5.68
        vcp = (
            1
            - (
                0.419869
                - 0.505839 * alpha
                - 1.56436 * alpha ** 3
                - 9481.7 * alpha ** 14
            )
        ) ** -8  # Eq 5.69
        sgp = (
            0.843593
            - 0.128624 * alpha
            - 3.36159 * alpha ** 3
            - 13749.5 * alpha ** 12
        )  # Eq 5.70
        return tb, tcp, pcp, vcp, sgp

    tb, Mp, tcp, pcp, vcp, sgp = Twu_tb(mw, sg, damp)
    #print(tb, Mp, tcp, pcp, vcp, sgp)
    tc = Twu_tc(tb, sgp, sg)

    vc = Twu_vc(tb, tcp, sg, sgp)
    pc = Twu_pc(tb, sgp, sg, pcp, tc, tcp, vc, vcp)
    #print(tc, pc, vc)
    if metric:
        tb = tb / 1.8  # Rankine -> Kelvin
        tc = tc / 1.8  # Rankine -> Kelvin
        pc = pc * PSI_TO_BAR  # psia -> barsa
    return (sg, tb, tc, pc, vc)

def sg_evolved_gas(
    p: float, degf: float, rsb: float, api: float, sg_sp: float
) -> float:
    """ Returns estimated specific gravity of gas evolved from oil insitu due to depressurization below Pb
        uses McCain & Hill Correlation (1995, SPE 30773)

        p: Pressure (psia)
        degf: Temperature (deg F)
        rsb: Oil solution GOR at Pb (scf/stb)
        api: Stock tank oil density (API)
        sg_sp: Specific gravity of separator gas (relative to air)
    """

    if (
        p > 314.7
    ):  # Two different sets from original 1995 paper (not reflected in Correlations book)
        a = [
            0,
            -208.0797,
            22885,
            -0.000063641,
            3.38346,
            -0.000992,
            -0.000081147,
            -0.001956,
            1.081956,
            0.394035,
        ]
    else:
        a = [
            0,
            -214.0887,
            9971,
            -0.001303,
            3.12715,
            -0.001495,
            -0.000085243,
            -0.003667,
            1.47156,
            0.714002,
        ]
    one_on_sgr = (
        a[1] / p
        + a[2] / p ** 2
        + a[3] * p
        + a[4] / degf ** 0.5
        + a[5] * degf
        + a[6] * rsb
        + a[7] * api
        + a[8] / sg_sp
        + a[9] * sg_sp ** 2
    )  # Eq 3.25
    return max(1 / one_on_sgr, sg_sp)


def sg_st_gas(
    psp: float, rsp: float, api: float, sg_sp: float, degf_sp: float
) -> float:
    """ Estimates specific gravity of gas evolving from stock tank
        from oil API and separator gas properties & conditions
        Returns sg_st (Stock Tank gas SG relative to air).
        Correlation reproduced from Valko McCain 2003 paper Eq 4-2

        psp: Separator pressure (psia)
        rsp: Separator GOR (separator scf / stb)
        api: Stock tank oil density (API)
        sg_sp: Separator gas specific gravity relative to air
        degf_sp: Separator temperature (deg f)
    """
    var = [np.log(psp), np.log(rsp), api, sg_sp, degf_sp]
    C = [
        [-17.275, -0.3354, 3.705, -155.52, 2.085],
        [7.9597, -0.3346, -0.4273, 629.61, -7.097e-2],
        [-1.1013, 0.1956, 1.818e-2, -957.38, 9.859e-4],
        [2.7735e-2, -3.4374e-2, -3.459e-4, 647.57, -6.312e-6],
        [3.2287e-3, 2.08e-3, 2.505e-6, -163.26, 1.4e-8],
    ]
    Zn = [sum([C[i][n] * var[n] ** i for i in range(5)]) for n in range(5)]
    Z = sum(Zn)
    sg_st = (
        1.219 + 0.198 * Z + 0.0845 * Z ** 2 + 0.03 * Z ** 3 + 0.003 * Z ** 4
    )
    return sg_st


def sgg_wt_avg(sg_sp: float, rsp: float, sg_st: float, rst: float) -> float:
    """ Calculates weighted average specific gravity of surface gas (sg_g)
        from separator and stock tank gas properties
        Returns sg_g (Weighted average surface gas SG relative to air).
        From McCain Correlations book, Eq 3.4

        sg_sp: Separator gas specific gravity relative to air
        rsp: Separator GOR (separator scf / stb)
        sg_st: Stock tank gas specific gravity relative to air
        rst: Stock tank producing gas-oil ratio (scf/stb)
    """
    sg_g = (sg_sp * rsp + sg_st * rst) / (rsp + rst)
    return sg_g


def oil_rs_st(psp: float, degf_sp: float, api: float, metric: bool = False) -> float:
    """ Estimates incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (scf/stb)

        Rsb = Rsp + Rst (Solution GOR at bubble point = Separator GOR + Stock Tank GOR).
        In absence of separator properties, a simple linear relationship with Rsp could be used instead;
          rs_st = 0.1618 * Separator GOR (Adapted from Eq 3-4 in Valko McCain 2003 paper)
        Correlation reproduced from Valko McCain 2003 paper Eq 3-2

        psp: Separator pressure (psia | barsa)
        degf_sp: Separator temperature (deg f | deg C)
        api: Stock tank oil density (API)
        metric: If True, input/output in Eclipse METRIC units. Defaults to False (FIELD)
    """
    if metric:
        psp = psp * BAR_TO_PSI
        degf_sp = degc_to_degf(degf_sp)

    var = [np.log(psp), np.log(degf_sp), api]
    C = [[-8.005, 1.224, -1.587], [2.7, -0.5, 0.0441], [-0.161, 0, -2.29e-5]]
    Zn = [sum([C[i][n] * var[n] ** i for i in range(3)]) for n in range(3)]
    Z = sum(Zn)
    return max(0, 3.955 + 0.83 * Z - 0.024 * Z ** 2 + 0.075 * Z ** 3)

def oil_pbub(
    api: float,
    degf: float,
    rsb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
) -> float:
    """ Returns bubble point pressure (psia | barsa) calculated with different correlations

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        rsb: Oil solution gas volume at Pbub (scf/stb | sm3/sm3)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997)
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        degf = degc_to_degf(degf)
        rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod = validate_methods(["pbmethod"], [pbmethod])

    if pbmethod.name == "STAN":
        if rsb * api * sg_g * degf == 0:
            raise ValueError(
                "Need valid values for rs, api, sg_g and degf for Standing Pb calculation"
            )
    else:
        if rsb * api * sg_sp * degf == 0:
            raise ValueError(
                "Need valid values for rsb, api, sg_sp and degf for Velarde or Valko McCain Pb calculation"
            )

    def pbub_standing(
        api, degf, sg_g, rsb, sg_sp
    ) -> float:  # 1.63 in 'Oil & Gas Properties & Correlations' - http://dx.doi.org/10.1016/B978-0-12-803437-8.00001-4
        a = 0.00091 * degf - 0.0125 * api
        return (
            18.2 * ((rsb / sg_g) ** 0.83 * 10 ** a - 1.4) + psc
        )  # Adding 14.7 as I suspect this is in psig

    def pbub_valko_mccain(api, degf, sg_g, rsb, sg_sp) -> float:
        extrap = False
        rsb2 = rsb
        if rsb <= 1: # Extrapolate to 14.696 psia at zero GOR
            extrap = True
            rsb = 1
        var = [np.log(rsb), api, sg_sp, degf]
        C = [
            [-5.48, 1.27, 4.51, -0.7835],
            [-0.0378, -0.0449, -10.84, 6.23e-3],
            [0.281, 4.36e-4, 8.39, -1.22e-5],
            [-0.0206, -4.76e-6, -2.34, 1.03e-8],
        ]
        Zn = [sum([C[i][n] * var[n] ** i for i in range(4)]) for n in range(4)] # Eq 2-1
        Z = sum(Zn)
        lnpb = 7.475 + 0.713 * Z + 0.0075 * Z ** 2 
        pb = np.exp(lnpb)
        
        if extrap:
            slope = (pb - 14.696)/(1 - 0)
            intercept = pb - 1*slope
            pb = slope * rsb2 + intercept

        return pb

    def pbub_velarde(api, degf, sg_g, rsb, sg_sp) -> float:
        x = 0.013098 * degf ** 0.282372 - 8.2e-6 * api ** 2.176124
        
        rsb_lim = (0.740152 / (sg_sp ** -0.161488 * 10 ** x))**(1/0.081465) # If Rsb < than this value, then the term inside the pbp brackets of pbp goes negative and causes imaginary numbers when raised to a power
        if rsb < rsb_lim+1e-6:
            rsb = rsb_lim+1e-6
            
        pbp = (
            1091.47
            * (rsb ** 0.081465 * sg_sp ** -0.161488 * 10 ** x - 0.740152)
            ** 5.354891 # psig
        )
        return pbp+psc # psia

    fn_dic = {
        "STAN": pbub_standing,
        "VALMC": pbub_valko_mccain,
        "VELAR": pbub_velarde,
    }

    result = fn_dic[pbmethod.name](
        api=api, degf=degf, sg_g=sg_g, rsb=rsb, sg_sp=sg_sp
    )
    if metric:
        return result * PSI_TO_BAR  # psia -> barsa
    return result

def oil_rs_bub(
    api: float,
    degf: float,
    pb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
    metric: bool = False,
) -> float:
    """ Returns Solution GOR (scf/stb | sm3/sm3) at bubble point pressure.
        Uses the inverse of the Bubble point pressure correlations, with the same method families
        Note: At low pressures, the VALMC method will fail (generally when Rsb < 10 scf/stb).
              The VALMC method will revert to the STAN method in these cases

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point Pressure (psia | barsa)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
        rsmethod: A string or pb_method Enum class that specifies one of following calculation choices. Note that VASBG is not available as it requires Pb as an input;
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
                   STAN: Standing Correlation (1947), using form from https://www.sciencedirect.com/science/article/pii/B9780128034378000014
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        degf = degc_to_degf(degf)
        pb = pb * BAR_TO_PSI

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod = 'VALMC' # Pbub calculations only needed with 'VALMC' rsmethod
    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )


    def rsbub_standing(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Standing')
        a = 0.00091 * degf - 0.0125 * api  # Eq 1.64
        return sg_g * (((pb - psc) / 18.2 + 1.4) / 10 ** a) ** (
            1 / 0.83
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def rsbub_valko_mccain(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Valko McCain')
        # Solve via iteration. First guess using Velarde Rsb, then simple Newton Iterations
        old_rsb = rsbub_velarde(api, degf, pb, sg_g, sg_sp)
        standing = False

        old_pbcalc = oil_pbub(degf=degf, api=api, sg_sp=sg_sp, rsb=old_rsb, pbmethod=pbmethod)
        old_err = old_pbcalc - pb
        new_rsb = old_rsb * pb / old_pbcalc
        i = 0
        new_err = 1000
        while abs(new_err) > 1e-5:
            i += 1
            new_pbcalc = oil_pbub(degf=degf,api=api,sg_sp=sg_sp,rsb=new_rsb,pbmethod=pbmethod)
            new_err = new_pbcalc - pb
            
            error_slope = (new_rsb - old_rsb) / (new_err - old_err)
            intcpt = new_rsb - error_slope * new_err
            
            old_err, old_rsb = new_err, new_rsb    
            
            new_rsb = intcpt
                
            #print('perr, pb_calc, rs_calc', new_err, new_pbcalc, new_rsb)
            if (i > 100):
                import warnings
                warnings.warn("oil_rs_bub: Valko-McCain did not converge after 100 iterations", RuntimeWarning)
                return new_rsb
        return new_rsb

    def rsbub_velarde(api, degf, pb, sg_g, sg_sp) -> float:
        x = 0.013098 * degf ** 0.282372 - (8.2e-6 * api ** 2.176124) # Eq 14 
        
        # Note that the Velarde approach predicts increasing Pbub with rs falling below 1 scf/stb, so below this we will linearly extrapolate to zero instead
        p_1scfstb = 1091.47*(sg_sp**-0.161488 * 10**x - 0.740152)**5.354891 # Eq 13
        psig = pb - 14.696        
        
        if psig >= p_1scfstb:
            rsb = (-10**(-x)* sg_sp**(0.161488) *(-(psig/1091.47)**(1/5.354891) - 0.740152))**(1/0.081465) # Rearranged Eq 13
        else:
            slope = (1-0)/(p_1scfstb - 0)
            intercept = 1 - slope * p_1scfstb
            rsb = slope * psig + intercept
        return rsb

    fn_dic = {
        "STAN": rsbub_standing,
        "VALMC": rsbub_valko_mccain,
        "VELAR": rsbub_velarde,
    }
    
    rsbub = fn_dic[rsmethod.name](
        api=api, degf=degf, pb=pb, sg_g=sg_g, sg_sp=sg_sp
    )
    if np.isnan(rsbub):
        import warnings
        warnings.warn("oil_rs_bub: correlation returned NaN, returning 0", RuntimeWarning)
        return 0
    if metric:
        return rsbub * SCF_PER_STB_TO_SM3_PER_SM3  # scf/stb -> sm3/sm3
    return rsbub

def oil_rs(
    api: float,
    degf: float,
    sg_sp: float,
    p: float,
    pb: float = 0,
    rsb: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
) -> float:
    """ Returns solution gas oil ratio (scf/stb | sm3/sm3) calculated from different correlations. Either pb, rsb or both need to be specified.
        If one is missing, the other will be calculated from correlation.

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        sg_sp: SG of separator gas
        p: Pressure of oil (psia | barsa)
        pb: Original bubble point pressure of oil (psia | barsa)
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb | sm3/sm3) <-- Required for Velarde, Blasingame & McCain
        rsmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
                   STAN: Standing Correlation (1947), using form from https://www.sciencedirect.com/science/article/pii/B9780128034378000014
                   VALMC: Valko-McCain Correlation (2003)
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        degf = degc_to_degf(degf)
        p = p * BAR_TO_PSI
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )

    sg_g, sg_sp = check_sgs(sg_g=0, sg_sp=sg_sp)

    if pb <= 0:  # Calculate Pb
        pb = oil_pbub(
            api=api, degf=degf, rsb=rsb, sg_sp=sg_sp, pbmethod=pbmethod
        )
    if rsb <= 0:  # Calculate rsb
        rsb = oil_rs_bub(
            api=api,
            degf=degf,
            pb=pb,
            sg_sp=sg_sp,
            rsmethod=rsmethod,
        )
        rsb = get_real_part(rsb)
        
    #print(rsb)
    
    if p >= pb:
        if metric:
            return rsb * SCF_PER_STB_TO_SM3_PER_SM3
        return rsb

    def Rs_velarde(
        api, degf, sg_sp, p, pb, rsb
    ):  # Velarde, Blasingame & McCain (1997)
        # Equations 3.8a - 3.8f
        # Estimates Rs of depleting oil from separator oil observations
        pb = max(psc, pb)
        p = max(psc, p)

        if sg_sp * api * rsb == 0:
            raise ValueError(
                "Missing one of the required inputs: sg_sp, api, rsb, for the Velarde, Blasingame & McCain Rs calculation"
            )
        A = [9.73e-7, 1.672608, 0.929870, 0.247235, 1.056052]
        B = [0.022339, -1.004750, 0.337711, 0.132795, 0.302065]
        C = [0.725167, -1.485480, -0.164741, -0.091330, 0.047094]

        xs = [A, B, C]
        a = [
            x[0]
            * sg_sp ** x[1]
            * api ** x[2]
            * degf ** x[3]
            * (pb - psc) ** x[4]
            for x in xs
        ]
        pr = (p - psc) / (pb - psc)
        rsr = a[0] * pr ** a[1] + (1 - a[0]) * pr ** a[2] # Eq 3 from Velarde & Blasingame
        rs = rsb * rsr
        return get_real_part(rs)

    def rs_standing(api, degf, sg_sp, p, pb, rsb):
        a = 0.00091 * degf - 0.0125 * api  # Eq 1.64
        return sg_g * (((p - psc) / 18.2 + 1.4) / 10 ** a) ** (
            1 / 0.83
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def rs_valko_mccain(api, degf, sg_sp, p, pb, rsb):
        rsb_valko = oil_rs_bub(api, degf, pb, sg_g, sg_sp, rsmethod = 'VALMC') # Rsb from Valko-McCain approach
        rs_scaler = rsb / rsb_valko # Scalar to adjust calculated Valko McCain Rs back to be in line with the supplied Rsb
        return rs_scaler * oil_rs_bub(api, degf, p, sg_g, sg_sp, rsmethod = 'VALMC')
            
                
    #def Rs_vasquezbegs(api, degf, sg_sp, p, pb, rsb):
    #    sg_gs = sg_sp * (
    #        1 + 5.912e-5 * api * degf_sep * np.log10(p_sep / 114.7)
    #    )  # Gas sg normalized to 100 psig separator conditions
    #    if api <= 30:
    #        return (
    #            0.0362
    #            * sg_gs
    #            * p ** 1.0937
    #            * np.exp(25.7240 * (api / (degf + degF2R)))
    #        )
    #    else:
    #        return (
    #            0.0178
    #            * sg_gs
    #            * p ** 1.1870
    #            * np.exp(23.9310 * (api / (degf + degF2R)))
    #        )

    fn_dic = {
        "VELAR": Rs_velarde,
        "STAN": rs_standing,
        #"VASBG": Rs_vasquezbegs,
        "VALMC": rs_valko_mccain,
    }

    result = fn_dic[rsmethod.name](
        api=api, degf=degf, sg_sp=sg_sp, p=p, pb=pb, rsb=rsb
    )
    if metric:
        return result * SCF_PER_STB_TO_SM3_PER_SM3  # scf/stb -> sm3/sm3
    return result

def check_sgs(
    sg_g: float,
    sg_sp: float,
    rst: float = 5,
    rsp: float = 1000,
    sg_st: float = 1.15,
) -> Tuple:
    """ Function used to impute sg_g or sg_sp when one or the other is zero
        sg_g: The weighted average surface-gas specific gravity (sep gas + gas evolved from liquid after separation)
        sg_sp: Separator specific gas gravity

        Optional
            rst: Post separation GOR (scf/stb)
            rsp: Separator GOR (scf/stb)
            sg_st: Gas sg evolved from post separartor liquid (rel air)

    """
    
    
    if sg_g > 0 and sg_sp > 0:
        return sg_g, sg_sp
    if sg_g <= 0 and sg_sp > 0:  # Estimate sg_g from sg_sp
        sg_g = (sg_sp * rsp + sg_st * rst) / (rsp + rst)
    if sg_g > 0 and sg_sp <= 0:  # Estimate sg_sp from sg_g
        sg_sp = (sg_g * (rsp + rst) - (sg_st * rst)) / rsp
    if sg_g < sg_sp:
        sg_sp = sg_g
    return (sg_g, sg_sp)

def oil_co(
    p: float,
    api: float,
    degf: float,
    sg_sp: float = 0,
    sg_g: float = 0,
    pb: float = 0,
    rsb: float = 0,
    pi: float = 0,
    comethod: co_method = co_method.EXPLT,
    zmethod: z_method = z_method.DAK,
    rsmethod: rs_method = rs_method.VELAR,
    cmethod: c_method = c_method.PMC,
    denomethod: deno_method = deno_method.SWMH,
    bomethod: bo_method = bo_method.MCAIN,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
):
    """ Returns oil compressibility (1/psi | 1/bar) calculated with Co = -1/Bo *[dBodp - Bg*dRsdp]
        using numerically derived values and their derivatives

        p: Reservoir pressure (psia | barsa)
        api: Stock tank oil density (deg API)
        sg_sp: Separator Gas specific Gravity (relative to air). If not defined, will use sg_g instead
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not defined, will use sg_sp instead
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point pressure (psia | barsa). If not provided, will attempt to calculate with Valko-McCain Pb Correlation
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb | sm3/sm3)
        comethod: A string or co_method Enum class that specifies calculation method for compressibility (currently only one option)
        zmethod: A string or z_method Enum class that specifies calculation method for gas Z-factor
        rsmethod: A string or rs_method Enum class that specifies calculation method for GOR
        cmethod: A string or c_method Enum class that specifies calculation method for gas critical properties
        denomethod: A string or deno_method Enum class that specifies calculation method for live oil density
        bomethod: A string or bo_method Enum class that specifies calculation method for oil FVF
        pbmethod: A string or pb_method Enum class that specifies calculation method for bubble point pressure
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3, 1/bar). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        if pi > 0:
            pi = pi * BAR_TO_PSI

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    (
        zmethod,
        rsmethod,
        cmethod,
        denomethod,
        bomethod,
        pbmethod,
        comethod,
    ) = validate_methods(
        [
            "zmethod",
            "rsmethod",
            "cmethod",
            "denomethod",
            "bomethod",
            "pbmethod",
            "comethod",
        ],
        [zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod, comethod],
    )

    if pb <= 0:  # Calculate Pb
        pb = oil_pbub(
            api=api, degf=degf, rsb=rsb, sg_sp=sg_sp, pbmethod=pbmethod
        )
    if rsb <= 0:  # Calculate rsb
        rsb = oil_rs_bub(
            api=api,
            degf=degf,
            pb=pb,
            sg_sp=sg_sp,
            rsmethod=rsmethod,
        )

    def Co_explicit(
        p=p,
        api=api,
        sg_sp=sg_sp,
        sg_g=sg_g,
        degf=degf,
        pb=pb,
        rsb=rsb,
        zmethod=zmethod,
        rsmethod=rsmethod,
        cmethod=cmethod,
        denomethod=denomethod,
        bomethod=bomethod,
    ):  # Explicit - Calculate with numerical derivatives
        # co = -1/bo*(dbodp - bg*drsdp/CUFTperBBL)
        def calc_dbodp(p):
            if p > pb:
                rs = rsb
            else:
                rs = oil_rs(
                    api=api,
                    degf=degf,
                    sg_sp=sg_sp,
                    p=p,
                    pb=pb,
                    rsb=rsb,
                    rsmethod=rsmethod,
                    pbmethod=pbmethod,
                )
            sg_o = oil_sg(api)
            bo = oil_bo(
                p=p,
                pb=pb,
                degf=degf,
                rs=rs,
                rsb=rsb,
                sg_sp=sg_sp,
                sg_g=sg_g,
                sg_o=sg_o,
                bomethod=bomethod,
                denomethod=denomethod,
            )
            return bo

        def calc_drsdp(p):
            return oil_rs(
                api=api,
                degf=degf,
                sg_sp=sg_sp,
                p=p,
                pb=pb,
                rsb=rsb,
                rsmethod=rsmethod,
                pbmethod=pbmethod,
            )

        dp = max(0.5, p * 0.001)  # Relative step size for numerical derivative
        if p > pb - dp and p < pb + dp:
            # Near bubble point, use one-sided derivative to avoid crossing Pb
            dbodp = calc_dbodp(p) - calc_dbodp(p - dp)
            drsdp = calc_drsdp(p) - calc_drsdp(p - dp)
        elif p > dp:
            dbodp = calc_dbodp(p + dp) - calc_dbodp(p - dp)
            drsdp = calc_drsdp(p + dp) - calc_drsdp(p - dp)
        else:
            dbodp = calc_dbodp(p + dp) - calc_dbodp(p)
            drsdp = calc_drsdp(p + dp) - calc_drsdp(p)

        if p > pb:
            drsdp = 0

        bo = calc_dbodp(p)
        bg = (
            gas.gas_bg(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
            / CUFTperBBL
        )  # rb/scf
        return -1 / bo * (dbodp - bg * drsdp)

    fn_dic = {"EXPLT": Co_explicit}

    result = fn_dic[comethod.name](
        p=p,
        api=api,
        sg_sp=sg_sp,
        sg_g=sg_g,
        degf=degf,
        pb=pb,
        rsb=rsb,
        zmethod=zmethod,
        rsmethod=rsmethod,
        cmethod=cmethod,
        denomethod=denomethod,
        bomethod=bomethod,
    )
    if metric:
        return result * INVPSI_TO_INVBAR  # 1/psi -> 1/bar
    return result

def oil_deno(
    p: float,
    degf: float,
    rs: float,
    rsb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    pb: float = 1e6,
    sg_o: float = 0,
    api: float = 0,
    denomethod: deno_method = deno_method.SWMH,
    metric: bool = False,
) -> float:
    """ Returns live oil density (lb/cuft | kg/m3) calculated with different correlations

        p: Pressure (psia | barsa)
        pb: Bubble point pressure (psia | barsa). Defaults to 1E6, and not used for densities below Pb. A valid value is required for density calculations above Pb
        degf: Reservoir Temperature (deg F | deg C)
        rs: Oil solution gas volume (scf/stb | sm3/sm3)
        rsb: Oil solution gas volume at bubble point pressure (scf/stb | sm3/sm3)
        sg_g: Weighted average specific gravity of surface gas (relative to air).
        sg_sp: Separator gas specific gravity (relative to air). If not known, an alternate nethod to estimate pseudo liquid density of surface gas will be used
        sg_o: Stock tank oil specific gravity (SG relative to water). If undefined will calculate from api
        api: Stock tank oil density (deg API). If undefined will calculate from sg_o. If both defined api value will prevail
        denomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   SWMH: Standing, Witte, McCain-Hill (1995) - Default
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3, kg/m3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        rs = rs * SM3_PER_SM3_TO_SCF_PER_STB
        rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        if pb != 1e6:
            pb = pb * BAR_TO_PSI

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    denomethod = validate_methods(["denomethod"], [denomethod])

    if sg_g == 0 and sg_sp == 0:
        raise ValueError(
            "Must define at least one of sg_g and sg_sp for density calculation"
        )

    # Density at or below initial bubble point pressure
    def Deno_standing_white_mccainhill(
        p: float,
        degf: float,
        rs: float,
        rsb: float,
        sg_g: float,
        sg_sp: float,
        pb: float,
        sg_o: float,
        api: float,
    ) -> float:  # (1995), Eq 3.18a - 3.18g

        if sg_sp > 0:
            a = np.array(
                [-49.8930, 85.0149, -3.70373, 0.0479818, 2.98914, -0.0356888]
            )
            rho_po = max(52.8 - 0.01 * rs, 20)  # First estimate
            err = 1
            i = 0
            while err > 1e-8:
                i += 1
                rhoa = (
                    a[0]
                    + a[1] * sg_sp
                    + a[2] * sg_sp * rho_po
                    + a[3] * sg_sp * rho_po ** 2
                    + a[4] * rho_po
                    + a[5] * rho_po ** 2
                )  # Eq 3.18c
                new_rho_po = (rs * sg_sp + 4600 * sg_o) / (
                    73.71 + rs * sg_sp / rhoa
                )  # pseudoliquid density, Eq 3.18b. Note equation in original paper uses sg_sp rather than sg_g as in book.
                err = abs(rho_po - new_rho_po)
                rho_po = new_rho_po
                if i > 100:
                    break
        else:
            rhoa = 38.52 * (10 ** (-0.00326 * api)) + (
                94.75 - 33.93 * np.log10(api)
            ) * np.log10(
                sg_g
            )  # Eq 3.17e using sg_g. Apparent liquid density of surface gases
            rho_po = (rs * sg_g + 4600 * sg_o) / (
                73.71 + rs * sg_g / rhoa
            )  # pseudoliquid density, Eq 3.18b

        drho_p = (
            (0.167 + 16.181 * 10 ** (-0.0425 * rho_po)) * p / 1000
            - 0.01 * (0.299 + 263 * 10 ** (-0.0603 * rho_po)) * (p / 1000) ** 2
        )  # Eq 3.19d

        rho_bs = rho_po + drho_p  # fake density used in calculations, Eq 3.19e
        drho_t = (
            (0.00302 + 1.505 * rho_bs ** -0.951) * (degf - tsc) ** 0.938
            - (0.0216 - 0.0233 * 10 ** (-0.0161 * rho_bs))
            * (degf - tsc) ** 0.475
        )  # Eq 3.19f
        rho_or = rho_bs - drho_t  # Eq 3.19g

        return rho_or

    def Deno_p_gt_pb(
        p: float,
        degf: float,
        rs: float,
        rsb: float,
        sg_g: float,
        sg_sp: float,
        pb: float,
        sg_o: float,
        api: float,
    ) -> float:
        rhorb = Deno_standing_white_mccainhill(
            p=pb,
            degf=degf,
            rs=rs,
            rsb=rsb,
            sg_g=sg_g,
            sg_sp=sg_sp,
            pb=pb,
            sg_o=sg_o,
            api=api,
        )

        # cofb calculation from default compressibility algorithm Eq 3.13
        C = [
            [3.011, -0.0835, 3.51, 0.327, -1.918, 2.52],
            [-2.6254, -0.259, -0.0289, -0.608, -0.642, -2.73],
            [0.497, 0.382, -0.0584, 0.0911, 0.154, 0.429],
        ]
        var = [
            np.log(api),
            np.log(sg_sp),
            np.log(pb),
            np.log(p / pb),
            np.log(rsb),
            np.log(degf),
        ]
        Zn = [sum([C[i][n] * var[n] ** i for i in range(3)]) for n in range(6)]
        Zp = sum(Zn)
        ln_cofb_p = 2.434 + 0.475 * Zp + 0.048 * Zp ** 2 - np.log(10 ** 6)
        cofb_p = np.exp(ln_cofb_p)

        return rhorb * np.exp(cofb_p * (p - pb))  # Eq 3.20

    fn_dic = {
        "SWMH": Deno_standing_white_mccainhill,
        "PGTPB": Deno_p_gt_pb,
    }  # Pressure greater than Pb

    if api == 0 and sg_o == 0:
        raise ValueError("Must supply either sg_o or api")

    if api == 0:  # Set api from sg_o
        api = 141.5 / sg_o - 131.5
    else:  # overwrite sg_o with api value
        sg_o = oil_sg(api)

    if (
        p > pb
    ):  # Use Eq 3.20, calculating oil density from density at Pb and compressibility factor
        result = fn_dic["PGTPB"](
            p=p,
            degf=degf,
            rs=rs,
            rsb=rsb,
            sg_g=sg_g,
            sg_sp=sg_sp,
            pb=pb,
            sg_o=sg_o,
            api=api,
        )
        if metric:
            return result * LBCUFT_TO_KGM3  # lb/cuft -> kg/m3
        return result

    result = fn_dic[denomethod.name](
        p=p,
        degf=degf,
        rs=rs,
        rsb=rsb,
        sg_g=sg_g,
        sg_sp=sg_sp,
        pb=pb,
        sg_o=sg_o,
        api=api,
    )
    if metric:
        return result * LBCUFT_TO_KGM3  # lb/cuft -> kg/m3
    return result

def oil_bo(
    p: float,
    pb: float,
    degf: float,
    rs: float,
    rsb: float,
    sg_o: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    bomethod: bo_method = bo_method.MCAIN,
    denomethod: deno_method = deno_method.SWMH,
    metric: bool = False,
) -> float:
    """ Returns oil formation volume factor (rb/stb | rm3/sm3) calculated with different correlations

        p: Pressure (psia | barsa)
        pb: Bubble point pressure (psia | barsa). Defaults to 1E6, and not used for densities below Pb. A valid value is required for density calculations above Pb
        degf: Reservoir Temperature (deg F | deg C)
        rs: Oil solution gas volume (scf/stb | sm3/sm3)
        rsb: Oil solution gas volume (scf/stb | sm3/sm3) at Pb
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not known, it will be estimated from sg_sp
        sg_sp: Separator gas specific gravity (relative to air).
        sg_o: Stock tank oil specific gravity (SG relative to water).
        bomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation
                   MCAIN: McCain approach, calculating from densities
        denomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   SWMH: Standing, Witte, McCain-Hill (1995) - Default
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        pb = pb * BAR_TO_PSI
        degf = degc_to_degf(degf)
        rs = rs * SM3_PER_SM3_TO_SCF_PER_STB
        rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)
    denomethod, bomethod = validate_methods(
        ["denomethod", "bomethod"], [denomethod, bomethod]
    )

    def Bo_standing(p, pb, degf, rs, rsb, sg_sp, sg_g, sg_o):
        Bob = (
            0.972
            + 1.47e-4 * (rs * (sg_g / sg_o) ** 0.5 + 1.25 * degf) ** 1.175
        )
        return Bob

    def Bo_mccain(p, pb, degf, rs, rsb, sg_sp, sg_g, sg_o):
        rhor = oil_deno(
            p=p,
            degf=degf,
            rs=rs,
            rsb=rsb,
            sg_g=sg_g,
            sg_sp=sg_sp,
            pb=pb,
            sg_o=sg_o,
            denomethod=denomethod,
        )
        return (sg_o * 62.372 + 0.01357 * rs * sg_g) / rhor  # Eq 3.21

    fn_dic = {"STAN": Bo_standing, "MCAIN": Bo_mccain}

    return fn_dic[bomethod.name](p, pb, degf, rs, rsb, sg_sp, sg_g, sg_o)


def oil_viso(p: float, api: float, degf: float, pb: float, rs: float, metric: bool = False) -> float:
    """ Returns Oil Viscosity (cP) with Beggs-Robinson (1975) correlation at saturated pressures
        and Petrosky-Farshad (1995) at undersaturated pressures

        p: Pressure (psia | barsa)
        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point Pressure (psia | barsa)
        rs: Solution GOR (scf/stb | sm3/sm3)
        metric: If True, input in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        pb = pb * BAR_TO_PSI
        rs = rs * SM3_PER_SM3_TO_SCF_PER_STB

    def uo_br(p, api, degf, pb, rs):
        Z = 3.0324 - 0.02023 * api
        y = 10 ** Z
        X = y * degf ** -1.163
        A = 10.715 * (rs + 100) ** -0.515
        B = 5.44 * (rs + 150) ** -0.338

        uod = 10 ** X - 1
        uor = A * uod ** B  # Eq 3.23c
        return uor

    def uo_pf(p, api, degf, pb, rs):
        uob = uo_br(pb, api, degf, pb, rs)
        loguob = np.log(uob)
        A = (
            -1.0146
            + 1.3322 * loguob
            - 0.4876 * loguob ** 2
            - 1.15036 * loguob ** 3
        )  # Eq 3.24b
        uor = uob + 1.3449e-3 * (p - pb) * 10 ** A  # Eq 3.24a
        return uor

    if p <= pb:
        return uo_br(p, api, degf, pb, rs)
    else:
        return uo_pf(p, api, degf, pb, rs)


class OilPVT:
    """ Oil PVT wrapper that stores oil characterization and method choices.
        Exposes rs(), bo(), density(), viscosity() methods delegating to oil_rs, oil_bo, oil_deno, oil_viso.

        api: Stock tank oil density (deg API)
        sg_sp: Separator gas specific gravity (relative to air)
        pb: Bubble point pressure (psia | barsa)
        rsb: Solution GOR at Pb (scf/stb | sm3/sm3)
        sg_g: Weighted average specific gravity of surface gas (relative to air). Estimated from sg_sp if not provided
        rsmethod: Method for Rs calculation. Defaults to 'VELAR'
        pbmethod: Method for Pb calculation. Defaults to 'VALMC'
        bomethod: Method for Bo calculation. Defaults to 'MCAIN'
        metric: If True, constructor inputs (pb, rsb) and method inputs/outputs use Eclipse METRIC units. Defaults to False (FIELD)
    """
    def __init__(self, api, sg_sp, pb, rsb, sg_g=0,
                 rsmethod='VELAR', pbmethod='VALMC', bomethod='MCAIN',
                 metric=False):
        self.api = api
        self.sg_sp = sg_sp
        self.metric = metric
        if metric:
            self.pb = pb * BAR_TO_PSI
            self.rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        else:
            self.pb = pb
            self.rsb = rsb
        self.sg_o = oil_sg(api)
        self.sg_g, self.sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)
        self.rsmethod = validate_methods(["rsmethod"], [rsmethod])
        self.pbmethod = validate_methods(["pbmethod"], [pbmethod])
        self.bomethod = validate_methods(["bomethod"], [bomethod])

    def rs(self, p, degf):
        """ Returns solution GOR (scf/stb | sm3/sm3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        if self.metric:
            p = p * BAR_TO_PSI
            degf = degc_to_degf(degf)
        result = oil_rs(api=self.api, degf=degf, sg_sp=self.sg_sp, p=p,
                      pb=self.pb, rsb=self.rsb, rsmethod=self.rsmethod,
                      pbmethod=self.pbmethod)
        if self.metric:
            return result * SCF_PER_STB_TO_SM3_PER_SM3
        return result

    def bo(self, p, degf, rs=None):
        """ Returns oil FVF (rb/stb | rm3/sm3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        if self.metric:
            p_field = p * BAR_TO_PSI
            degf_field = degc_to_degf(degf)
        else:
            p_field = p
            degf_field = degf
        if rs is None:
            rs_field = oil_rs(api=self.api, degf=degf_field, sg_sp=self.sg_sp, p=p_field,
                              pb=self.pb, rsb=self.rsb, rsmethod=self.rsmethod,
                              pbmethod=self.pbmethod)
        else:
            rs_field = rs * SM3_PER_SM3_TO_SCF_PER_STB if self.metric else rs
        return oil_bo(p=p_field, pb=self.pb, degf=degf_field, rs=rs_field, rsb=self.rsb,
                      sg_o=self.sg_o, sg_g=self.sg_g, sg_sp=self.sg_sp,
                      bomethod=self.bomethod)

    def density(self, p, degf, rs=None):
        """ Returns live oil density (lb/cuft | kg/m3) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        if self.metric:
            p_field = p * BAR_TO_PSI
            degf_field = degc_to_degf(degf)
        else:
            p_field = p
            degf_field = degf
        if rs is None:
            rs_field = oil_rs(api=self.api, degf=degf_field, sg_sp=self.sg_sp, p=p_field,
                              pb=self.pb, rsb=self.rsb, rsmethod=self.rsmethod,
                              pbmethod=self.pbmethod)
        else:
            rs_field = rs * SM3_PER_SM3_TO_SCF_PER_STB if self.metric else rs
        result = oil_deno(p=p_field, degf=degf_field, rs=rs_field, rsb=self.rsb,
                        sg_g=self.sg_g, sg_sp=self.sg_sp, pb=self.pb,
                        sg_o=self.sg_o)
        if self.metric:
            return result * LBCUFT_TO_KGM3
        return result

    def viscosity(self, p, degf, rs=None):
        """ Returns oil viscosity (cP) at pressure p (psia | barsa) and temperature degf (deg F | deg C) """
        if self.metric:
            p_field = p * BAR_TO_PSI
            degf_field = degc_to_degf(degf)
        else:
            p_field = p
            degf_field = degf
        if rs is None:
            rs_field = oil_rs(api=self.api, degf=degf_field, sg_sp=self.sg_sp, p=p_field,
                              pb=self.pb, rsb=self.rsb, rsmethod=self.rsmethod,
                              pbmethod=self.pbmethod)
        else:
            rs_field = rs * SM3_PER_SM3_TO_SCF_PER_STB if self.metric else rs
        return oil_viso(p=p_field, api=self.api, degf=degf_field, pb=self.pb, rs=rs_field)


def oil_harmonize_pb_rsb(
    pb: float = 0,
    rsb: float = 0,
    degf: float = 0,
    api: float = 0,
    sg_sp: float = 0,
    sg_g: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
    pbmethod: pb_method = pb_method.VELAR,
    metric: bool = False,
) -> Tuple:
    """Resolves consistent Pb, Rsb, and rsb_frac from user inputs.

    Given one or both of Pb and Rsb, returns values that are mutually consistent
    with the selected oil PVT correlations.

    - If only pb is specified (rsb=0): calculates rsb from pb
    - If only rsb is specified (pb=0): calculates pb from rsb
    - If both are specified: finds rsb_frac scaling factor that honors both values

    pb: Bubble point pressure (psia | barsa). Default 0 (unknown)
    rsb: Solution GOR at Pb (scf/stb | sm3/sm3). Default 0 (unknown)
    degf: Reservoir temperature (deg F | deg C)
    api: Stock tank oil density (deg API)
    sg_sp: Separator gas specific gravity
    sg_g: Weighted average surface gas specific gravity
    rsmethod: Rs calculation method. Default VELAR
    pbmethod: Pb calculation method. Default VELAR
    metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)

    Returns tuple of (pb, rsb, rsb_frac) where:
      - pb: Bubble point pressure (psia | barsa)
      - rsb: Solution GOR at bubble point (scf/stb | sm3/sm3)
      - rsb_frac: Scaling factor applied to Rs correlation to honor user Pb and Rsb
                  (1.0 if only one was specified)
    """
    if metric:
        degf = degc_to_degf(degf)
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    rsmethod, pbmethod = validate_methods(
        ["rsmethod", "pbmethod"], [rsmethod, pbmethod]
    )

    rsb_i = rsb
    pb_i = pb
    rsb_frac = 1.0

    # Calculate rsb from pb
    if rsb_i <= 0 and pb_i > 0:
        rsb = oil_rs_bub(
            degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
            pb=pb, rsmethod=rsmethod,
        )

    # Calculate pb from rsb
    if pb_i <= 0 and rsb_i > 0:
        pb = oil_pbub(
            degf=degf, api=api, sg_sp=sg_g, rsb=rsb, pbmethod=pbmethod
        )

    # Both defined by user — find rsb_frac that honors both
    if pb_i > 0 and rsb_i > 0:
        pbcalc = oil_pbub(
            degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
            rsb=rsb, pbmethod=pbmethod,
        )
        err = 100
        rsb_old = rsb
        i = 0
        while err > 0.0001:
            rsbnew = pb / pbcalc * rsb_old
            pbcalc = oil_pbub(
                degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
                rsb=rsbnew, pbmethod=pbmethod,
            )
            rsb_old = rsbnew
            err = abs(pb - pbcalc)
            i += 1
            if i > 100:
                raise RuntimeError(
                    "Could not solve Pb & Rsb for these combination of inputs"
                )
        rsb_frac = rsb_i / rsbnew

    if metric:
        return pb * PSI_TO_BAR, rsb * SCF_PER_STB_TO_SM3_PER_SM3, rsb_frac
    return pb, rsb, rsb_frac


def _resolve_pb_rsb(pb, rsb, degf, api, sg_sp, sg_g, pvto, pmax,
                     rsmethod, pbmethod):
    """Internal helper: resolves Pb/Rsb using oil_harmonize_pb_rsb plus PVTO extension.

    Returns (pb, rsb, rsb_frac, rsb_max, pb_i, rsb_i).
    """
    pb_i = pb
    rsb_i = rsb

    pb, rsb, rsb_frac = oil_harmonize_pb_rsb(
        pb=pb, rsb=rsb, degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
        rsmethod=rsmethod, pbmethod=pbmethod,
    )

    rsb_max = rsb

    if pvto and pmax > pb:
        rsb_max = oil_rs_bub(
            degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
            pb=pmax, rsmethod=rsmethod,
        )
        rs_at_pbi = oil_rs(
            api=api, degf=degf, sg_sp=sg_sp, p=pb, pb=pmax,
            rsb=rsb_max, rsmethod=rsmethod, pbmethod=pbmethod,
        )
        err = rs_at_pbi - rsb
        rsb_frac_new = rsb_frac
        i = 0
        while abs(err) > 0.0001:
            rsb_frac_new = rsb / rs_at_pbi * rsb_frac
            rs_at_pbi = oil_rs(
                api=api, degf=degf, sg_sp=sg_sp, p=pb, pb=pmax,
                rsb=rsb_max * rsb_frac_new,
                rsmethod=rsmethod, pbmethod=pbmethod,
            )
            rsb_frac = rsb_frac_new
            err = rs_at_pbi - rsb
            i += 1
            if i > 100:
                raise RuntimeError(
                    "Could not solve Pb & Rsb for these combination of inputs"
                )
        rsb_frac = rsb_frac_new

    return pb, rsb, rsb_frac, rsb_max, pb_i, rsb_i


def _build_bot_tables(pressures, pb, rsb, rsb_frac, rsb_max, sg_o, sg_g, sg_sp,
                      api, degf, pvto, wt, ch4_sat,
                      zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod):
    """Compute all PVT properties over the pressure array.

    Returns (rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
             usat_p, usat_bo, usat_uo) where usat_* are empty lists if not pvto.
    """
    if pvto:
        pb = max(pb, max(pressures))
        rsb = rsb_max * rsb_frac

    co, cg, rss, bos, uos, gfvf, visg, gz, bws, visws, denos = [
        [] for _ in range(11)
    ]

    for p in pressures:
        if p > pb:
            rss.append(rsb)
        else:
            rss.append(
                oil_rs(
                    api=api, degf=degf, sg_sp=sg_sp, p=p, pb=pb,
                    rsb=rsb / rsb_frac, rsmethod=rsmethod, pbmethod=pbmethod,
                ) * rsb_frac
            )

        bos.append(
            oil_bo(
                p=p, pb=pb, degf=degf, rs=rss[-1], rsb=rsb,
                sg_g=sg_g, sg_sp=sg_sp, sg_o=sg_o,
                denomethod=denomethod, bomethod=bomethod,
            )
        )
        denos.append(
            oil_deno(
                p=p, degf=degf, rs=rss[-1], rsb=rsb,
                sg_g=sg_g, sg_sp=sg_sp, pb=pb, sg_o=sg_o, api=api,
            )
        )
        uos.append(oil_viso(p=p, api=api, degf=degf, pb=pb, rs=rss[-1]))
        co.append(
            oil_co(
                p=p, api=api, sg_sp=sg_sp, sg_g=sg_g, degf=degf,
                pb=pb, rsb=rss[-1], zmethod=zmethod, rsmethod=rsmethod,
                cmethod=cmethod, denomethod=denomethod, bomethod=bomethod,
            )
        )

        gfvf.append(
            gas.gas_bg(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
            * 1000 / CUFTperBBL
        )
        gz.append(
            gas.gas_z(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        visg.append(
            gas.gas_ug(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        cg.append(gas.gas_cg(p=p, sg=sg_g, degf=degf, cmethod=cmethod))
        bw, _lden, visw, _cw, _rsw = brine.brine_props(
            p=p, degf=degf, wt=wt, ch4_sat=ch4_sat
        )
        bws.append(bw)
        visws.append(visw)

    # Undersaturated extension for PVTO
    usat_p, usat_bo, usat_uo = [], [], []
    if pvto:
        for i, p in enumerate(pressures):
            if i == 0:
                continue
            try:
                usat_p.append(pressures[i:])
                usat_bo.append(
                    [
                        oil_bo(
                            p=pusat, pb=p, degf=degf, rs=rss[i], rsb=rss[i],
                            sg_g=sg_g, sg_sp=sg_sp, sg_o=sg_o,
                            denomethod=denomethod, bomethod=bomethod,
                        )
                        for pusat in usat_p[-1]
                    ]
                )
                usat_uo.append(
                    [
                        oil_viso(p=pusat, api=api, degf=degf, pb=p, rs=rss[i])
                        for pusat in usat_p[-1]
                    ]
                )
            except (ValueError, IndexError, ZeroDivisionError):
                pass

    return (rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
            usat_p, usat_bo, usat_uo)


def _format_bot_results(pressures, rss, bos, denos, uos, co, gz, gfvf, cg,
                        visg, bws, visws, usat_p, usat_bo, usat_uo,
                        sg_o, sg_g, pi, degf, wt, ch4_sat, pb_i, rsb_i,
                        rsb_frac, pvto, export, zmethod, cmethod):
    """Assemble DataFrame, optionally export Eclipse files, and return results dict."""
    st_deno = sg_o * WDEN
    st_deng = gas.gas_den(
        p=psc, sg=sg_g, degf=tsc, zmethod=zmethod, cmethod=cmethod
    )
    bw, lden, visw, cw, _rsw = brine.brine_props(
        p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat
    )
    res_denw = lden * WDEN
    res_cw = cw

    df = pd.DataFrame()
    df["Pressure (psia)"] = pressures
    df["Rs (mscf/stb)"] = rss
    df["Rs (mscf/stb)"] = df["Rs (mscf/stb)"] / 1000
    df["Bo (rb/stb)"] = bos
    df["Deno (lb/cuft)"] = denos
    df["uo (cP)"] = uos
    df["Co (1/psi)"] = co
    df["Gas Z (v/v)"] = gz
    df["Bg (rb/mscf"] = gfvf
    df["Cg (1/psi)"] = cg
    df["ug (cP)"] = visg
    df["Bw (rb/stb)"] = bws
    df["uw (cP)"] = visws

    if export:
        df.to_excel("bot.xlsx", index=False, engine="openpyxl")
        pvdg = df[["Pressure (psia)", "Bg (rb/mscf", "ug (cP)"]]
        pvdg = pvdg.set_index("Pressure (psia)")
        headers = ["-- P (psia)", "Bg (rb/mscf", "ug (cP)"]
        fileout = "PVDG\n" + tabulate(pvdg, headers) + "\n/"
        with open("PVDG.INC", "w") as text_file:
            text_file.write(fileout)
        pvdo = df[["Pressure (psia)", "Bo (rb/stb)", "uo (cP)"]]
        pvdo = pvdo.set_index("Pressure (psia)")
        headers = ["-- P (psia)", "Bo (rb/stb)", "uo (cP)"]
        fileout = "PVDO\n" + tabulate(pvdo, headers) + "\n/"
        with open("PVDO.INC", "w") as text_file:
            text_file.write(fileout)

        if pvto:
            pvto_out = "PVTO\n"
            headers = [
                "-- Rs (mscf/stb)",
                "P (psia)",
                "Bo (rb/stb)",
                "uo (cP)",
                "",
            ]
            table = []
            for r, row in df.iterrows():
                table.append(
                    [
                        row["Rs (mscf/stb)"],
                        row["Pressure (psia)"],
                        row["Bo (rb/stb)"],
                        row["uo (cP)"],
                        "/",
                    ]
                )
                try:
                    if r > 0:
                        for e, entry in enumerate(usat_p[r - 1]):
                            if e == 0:
                                continue
                            table.append(
                                [
                                    " ",
                                    entry,
                                    usat_bo[r - 1][e],
                                    usat_uo[r - 1][e],
                                    " ",
                                ]
                            )
                except (IndexError, KeyError):
                    pass
            pvto_out += tabulate(table, headers)
            pvto_out += "\n/"
            with open("PVTO.INC", "w") as text_file:
                text_file.write(pvto_out)

    results = {
        "bot": df,
        "deno": st_deno,
        "deng": st_deng,
        "denw": res_denw,
        "cw": res_cw,
        "uw": visw,
        "pb": pb_i,
        "rsb": rsb_i,
        "rsb_scale": rsb_frac,
        "usat": [],
    }
    if pvto:
        results["usat"] = [usat_p, usat_bo, usat_uo]

    return results


def make_bot_og(*args, **kwargs):
    """Deprecated: Use simtools.make_bot_og() instead. This wrapper remains for backward compatibility."""
    from pyrestoolbox.simtools import make_bot_og as _make_bot_og
    return _make_bot_og(*args, **kwargs)
