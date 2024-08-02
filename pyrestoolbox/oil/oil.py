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
"""

import sys
from collections import Counter
import glob
from enum import Enum
import pkg_resources

import numpy as np
import numpy.typing as npt
import pandas as pd
from tabulate import tabulate

from pyrestoolbox.constants import R, psc, tsc, degF2R, tscr, MW_AIR, scf_per_mol, CUFTperBBL, WDEN
from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic
from pyrestoolbox.validate import validate_methods
import pyrestoolbox.gas as gas
import pyrestoolbox.brine as brine

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
) -> np.ndarray:
    """ Returns liquid rate for radial flow (stb/day) using Darcy pseudo steady state equation
        k: Effective Permeability to flow (mD)
        h: Net flow height (ft)
        Pr: Reservoir pressure (psia)
        pwf: BHFP (psia)
        r_w: Wellbore Radius (ft)
        r_ext: External Reservoir Radius (ft)
        uo: Liquid viscosity (cP)
        bo: Liquid Formation Volume Factor (rb/stb)
        S: Wellbore Skin (Dimensionless). Defaults to zero if not specified
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
        pb: Bubble point pressure (psia). Defaults to zero if not specified. Not used unless Vogel option is invoked
    """
    k, h, pr, pwf = (
        np.asarray(k),
        np.asarray(h),
        np.asarray(pr),
        np.asarray(pwf),
    )

    if pwf.size > 1:
        pb = np.array([max(pb, pwf[i]) for i in range(pwf.size)])

    J = (
        0.00708 * k * h / (uo * bo * (np.log(r_ext / r_w) + S - 0.75))
    )  # Productivity index
    if not vogel:
        qoil = J * (pr - pwf)
    else:
        qsat_max = J * pb / 1.8
        qusat = J * (pr - pb)
        qoil = (
            qsat_max * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2) + qusat
        )
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
) -> np.ndarray:
    """ Returns liquid rate for linear flow (stb/day) using Darcy steady state equation
        k: Permeability (mD)
        Pr: Reservoir pressure (psia)
        pwf: BHFP (psia)
        area: Net cross sectional area perpendicular to direction of flow (ft2).
        length: Length over which flow takes place (ft)
        uo: Liquid viscosity (cP)
        bo: Liquid Formation Volume Factor (rb/stb)
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
        pb: Bubble point pressure (psia). Defaults to zero if not specified. Not used unless Vogel option is invoked
    """
    k, area, pr, pwf = (
        np.asarray(k),
        np.asarray(area),
        np.asarray(pr),
        np.asarray(pwf),
    )

    J = (
        0.00708 * k * area / (2 * np.pi * uo * bo * length)
    )  # Productivity index
    if not vogel:
        qoil = J * (pr - pwf)
    else:
        qsat_max = J * pb / 1.8
        qusat = J * (pr - pb)
        qoil = (
            qsat_max * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2) + qusat
        )
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
    mw: float, ja: float = 0, sg: float = 0, damp: float = 1
) -> tuple:
    """ Returns tuple of sg, tb (degR), tc (DegR), pc (psia), Vc (ft3/lb-mol) using method from Twu (1984) correlations for petroleum liquids
        Modified with damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu)
        Returns sg, tb (R), tc (R), pc (psia), vc (ft3/lbmol)

        mw: Molecular weight of the liquid hydrocarbon (g/g.mol / lb/lb.mol)
        ja: jacoby Aromaticity Factor relationship. Varies between 0 (Paraffins) - 1 (Aromatic). Defaults to zero if undefined
        sg: Specific gravity of the liquid (fraction relative to water density). Will use jacoby method to estimate sg from mw if undefined.
        damp: damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu). Defaults to 1
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
                print("Check inputs. Twu algorithm did not converge", mw, ja, sg, damp)
                break
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


def oil_rs_st(psp: float, degf_sp: float, api: float) -> float:
    """ Estimates incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (scf/stb)

        Rsb = Rsp + Rst (Solution GOR at bubble point = Separator GOR + Stock Tank GOR).
        In absence of separator properties, a simple linear relationship with Rsp could be used instead;
          rs_st = 0.1618 * Separator GOR (Adapted from Eq 3-4 in Valko McCain 2003 paper)
        Correlation reproduced from Valko McCain 2003 paper Eq 3-2

        psp: Separator pressure (psia)
        degf_sp: Separator temperature (deg f)
        api: Stock tank oil density (API)
    """
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
) -> float:
    """ Returns bubble point pressure (psia) calculated with different correlations

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F)
        rsb: Oil solution gas volume at Pbub (scf/stb)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
    """
    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod = validate_methods(["pbmethod"], [pbmethod])

    if pbmethod.name == "STAN":
        if rsb * api * sg_g * degf == 0:
            print(
                "Need valid values for rs, api, sg_g for Standing Pb calculation"
            )
            print(
                "Need valid values for rs, api, sg_g and degf for Standing or Velarde Pb calculation"
            )
            sys.exit()
    else:
        if rsb * api * sg_sp * degf == 0:
            #print(rsb, api, sg_sp, degf)
            print(
                "Need valid values for rsb, api, sg_sp and degf for Velarde or Valko McCain Pb calculation"
            )
            sys.exit()

    def pbub_standing(
        api, degf, sg_g, rsb, sg_sp
    ) -> float:  # 1.63 in 'Oil & Gas Properties & Correlations' - http://dx.doi.org/10.1016/B978-0-12-803437-8.00001-4
        a = 0.00091 * degf - 0.0125 * api
        return (
            18.2 * ((rsb / sg_g) ** 0.83 * 10 ** a - 1.4) + psc
        )  # Adding 14.7 as I suspect this is in psig

    def pbub_valko_mccain(api, degf, sg_g, rsb, sg_sp) -> float:
        if rsb <= 10:
            return pbub_velarde(api, degf, sg_g, 0, sg_sp)
        var = [np.log(rsb), api, sg_sp, degf]
        C = [
            [-5.48, 1.27, 4.51, -0.7835],
            [-0.0378, -0.0449, -10.84, 6.23e-3],
            [0.281, 4.36e-4, 8.39, -1.22e-5],
            [-0.0206, -4.76e-6, -2.34, 1.03e-8],
        ]
        Zn = [sum([C[i][n] * var[n] ** i for i in range(4)]) for n in range(4)]
        Z = sum(Zn)
        lnpb = 7.475 + 0.713 * Z + 0.0075 * Z ** 2
        return np.exp(lnpb)

    def pbub_velarde(api, degf, sg_g, rsb, sg_sp) -> float:
        x = 0.013098 * degf ** 0.282372 - 8.2e-6 * api ** 2.176124
        pbp = (
            1091.47
            * (rsb ** 0.081465 * sg_sp ** -0.161488 * 10 ** x - 0.740152)
            ** 5.354891
        )
        return pbp

    fn_dic = {
        "STAN": pbub_standing,
        "VALMC": pbub_valko_mccain,
        "VELAR": pbub_velarde,
    }

    return fn_dic[pbmethod.name](
        api=api, degf=degf, sg_g=sg_g, rsb=rsb, sg_sp=sg_sp
    )

def oil_rs_bub(
    api: float,
    degf: float,
    pb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
) -> float:
    """ Returns Solution GOR (scf/stb) at bubble point pressure.
        Uses the inverse of the Bubble point pressure correlations, with the same method families
        Note: At low pressures, the VALMC method will fail (generally when Rsb < 10 scf/stb).
              The VALMC method will revert to the STAN method in these cases

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F)
        pb: Bubble point Pressure (psia)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
        rsmethod: A string or pb_method Enum class that specifies one of following calculation choices. Note that VASBG is not available as it requires Pb as an input;
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
                   STAN: Standing Correlation (1947), using form from https://www.sciencedirect.com/science/article/pii/B9780128034378000014
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
    """
    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod = 'VALMC' # Pbub calculations only needed with 'VALMC' rsmethod
    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )


    def rsbub_standing(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Standing')
        a = 0.00091 * degf - 0.0125 * api  # Eq 1.64
        return sg_sp * (((pb - psc) / 18.2 + 1.4) / 10 ** a) ** (
            1 / 0.83
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def rsbub_valko_mccain(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Valko McCain')
        # Solve via iteration. First guess using Velarde Rsb, then simple Newton Iterations
        old_rsb = rsbub_velarde(api, degf, pb, sg_g, sg_sp)
        standing = False
        if old_rsb < 10:
            return old_rsb
        old_pbcalc = oil_pbub(
            degf=degf, api=api, sg_sp=sg_sp, rsb=old_rsb, pbmethod=pbmethod
        )
        old_err = old_pbcalc - pb
        new_rsb = old_rsb * pb / old_pbcalc
        i = 0
        new_err = 1000
        while abs(new_err) > 1e-5:
            i += 1
            if new_rsb > 10:
                new_pbcalc = oil_pbub(
                    degf=degf,
                    api=api,
                    sg_sp=sg_sp,
                    rsb=new_rsb,
                    pbmethod=pbmethod,
                )
            else:
                new_pbcalc = oil_pbub(
                    degf=degf,
                    api=api,
                    sg_sp=sg_sp,
                    rsb=new_rsb,
                    pbmethod="STAN",
                )
            new_err = new_pbcalc - pb
            
            error_slope = (new_rsb - old_rsb) / (new_err - old_err)
            intcpt = new_rsb - error_slope * new_err
            
            old_err, old_rsb = new_err, new_rsb    
            
            new_rsb = intcpt
                
            #print('perr, pb_calc, rs_calc', new_err, new_pbcalc, new_rsb)
            if (i > 100):  # At low rsb VALMC will not converge, use Velarde instead
                return rsbub_velarde(api, degf, pb, sg_g, sg_sp)
        return new_rsb

    def rsbub_velarde(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Velarde')
        x = 0.013098 * degf ** 0.282372 - 8.2e-6 * api ** 2.176124
        rsb = (
            0.270811 * sg_sp ** (10093 / 62500) * pb ** 0.186745 * 10 ** (-x)
            + 92519 * sg_sp ** (10093 / 62500) * 2 ** (-x - 3) * 5 ** (-x - 6)
        ) ** (200000 / 16293)
        #print(api, degf, pb, sg_g, sg_sp, rsb)
        return max(rsb, 0)

    fn_dic = {
        "STAN": rsbub_standing,
        "VALMC": rsbub_valko_mccain,
        "VELAR": rsbub_velarde,
    }
    
    rsbub = fn_dic[pbmethod.name](
        api=api, degf=degf, pb=pb, sg_g=sg_g, sg_sp=sg_sp
    )
    if np.isnan(rsbub):
        return 0
    else:
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
) -> float:
    """ Returns solution gas oil ratio (scf/stb) calculated from different correlations. Either pb, rsb or both need to be specified.
        If one is missing, the other will be calculated from correlation.

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F)
        sg_sp: SG of separator gas
        p: Pressure of oil (psia)
        pb: Original bubble point pressure of oil (psia)
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb) <-- Required for Velarde, Blasingame & McCain
        rsmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
                   STAN: Standing Correlation (1947), using form from https://www.sciencedirect.com/science/article/pii/B9780128034378000014
                   VASBG: Vasquez & Beggs Correlation (1984)
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
    """
    #print(pbmethod, rsmethod)
    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )
    
    #print(sg_g, sg_sp, api, rsb)

    if pb <= 0:  # Calculate Pb
        pb = oil_pbub(
            api=api, degf=degf, rsb=rsb, sg_sp=sg_sp, pbmethod=pbmethod
        )
    if rsb <= 0:  # Calculate rsb
        #print('Calculating Rsb')
        #print(api, degf, pb, sg_sp)
        rsb = oil_rs_bub(
            api=api,
            degf=degf,
            pb=pb,
            sg_sp=sg_sp,
            rsmethod=rsmethod,
        )
        
    #print(rsb)
    
    if p > pb:
        return rsb

    def Rs_velarde(
        api, degf, sg_sp, p, pb, rsb
    ):  # Velarde, Blasingame & McCain (1997)
        # Equations 3.8a - 3.8f
        # Estimates Rs of depleting oil from separator oil observations
        pb = max(psc, pb)
        p = max(psc, p)

        if sg_sp * api * rsb == 0:
            print(
                "Missing one of the required inputs: sg_sp, api, rsb, for the Velarde, Blasingame & McCain Rs calculation"
            )
            sys.exit()
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
        rsr = a[0] * pr ** a[1] + (1 - a[0]) * pr ** a[2]
        rs = rsb * rsr
        return rs

    def rs_standing(api, degf, sg_sp, p, pb, rsb):
        a = 0.00091 * degf - 0.0125 * api  # Eq 1.64
        return sg_sp * (((p - psc) / 18.2 + 1.4) / 10 ** a) ** (
            1 / 0.83
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def Rs_vasquezbegs(api, degf, sg_sp, p, pb, rsb):
        sg_gs = sg_sp * (
            1 + 5.912e-5 * api * degf_sep * np.log10(p_sep / 114.7)
        )  # Gas sg normalized to 100 psig separator conditions
        if api <= 30:
            return (
                0.0362
                * sg_gs
                * p ** 1.0937
                * np.exp(25.7240 * (api / (degf + degF2R)))
            )
        else:
            return (
                0.0178
                * sg_gs
                * p ** 1.1870
                * np.exp(23.9310 * (api / (degf + degF2R)))
            )

    fn_dic = {
        "VELAR": Rs_velarde,
        "STAN": rs_standing,
        "VASBG": Rs_vasquezbegs,
    }

    return fn_dic[rsmethod.name](
        api=api, degf=degf, sg_sp=sg_sp, p=p, pb=pb, rsb=rsb
    )

def check_sgs(
    sg_g: float,
    sg_sp: float,
    rst: float = 5,
    rsp: float = 1000,
    sg_st: float = 1.15,
) -> tuple:
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
):
    """ Returns oil compressibility (1/psi) calculated with Co = -1/Bo *[dBodp - Bg*dRsdp]
        using numerically derived values and their derivatives

        p: Reservoir pressure (psia)
        api: Stock tank oil density (deg API)
        sg_sp: Separator Gas specific Gravity (relative to air). If not defined, will use sg_g instead
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not defined, will use sg_sp instead
        degf: Reservoir Temperature (deg F)
        pb: Bubble point pressure. If not provided, will attempt to calculate with Valko-McCain Pb Correlation
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb)
        comethod: A string or co_method Enum class that specifies calculation method for compressibility (currently only one option)
        zmethod: A string or z_method Enum class that specifies calculation method for gas Z-factor
        rsmethod: A string or rs_method Enum class that specifies calculation method for GOR
        cmethod: A string or c_method Enum class that specifies calculation method for gas critical properties
        denomethod: A string or deno_method Enum class that specifies calculation method for live oil density
        bomethod: A string or bo_method Enum class that specifies calculation method for oil FVF
        pbmethod: A string or pb_method Enum class that specifies calculation method for bubble point pressure
    """
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

        if p > 15.7:
            if p < pb - 0.5 or p > pb + 0.5:
                dbodp = calc_dbodp(p + 0.5) - calc_dbodp(p - 0.5)
                drsdp = calc_drsdp(p + 0.5) - calc_drsdp(p - 0.5)
            else:
                dbodp = calc_dbodp(p) - calc_dbodp(p - 1)
                drsdp = calc_drsdp(p) - calc_drsdp(p - 1)
        else:
            dbodp = calc_dbodp(p + 1) - calc_dbodp(p)
            drsdp = calc_drsdp(p + 1) - calc_drsdp(p)

        if p > pb:
            drsdp = 0

        bo = calc_dbodp(p)
        bg = (
            gas.gas_bg(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
            / CUFTperBBL
        )  # rb/scf
        return -1 / bo * (dbodp - bg * drsdp)

    fn_dic = {"EXPLT": Co_explicit}

    return fn_dic[comethod.name](
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
) -> float:
    """ Returns live oil density calculated with different correlations

        p: Pressure (psia)
        pb: Bubble point pressure (psia). Defaults to 1E6, and not used for densities below Pb. A valid value is required for density calculations above Pb
        degf: Reservoir Temperature (deg F)
        rs: Oil solution gas volume (scf/stb)
        rsb: Oil solution gas volume at bubble point pressure (scf/stb)
        sg_g: Weighted average specific gravity of surface gas (relative to air).
        sg_sp: Separator gas specific gravity (relative to air). If not known, an alternate nethod to estimate pseudo liquid density of surface gas will be used
        sg_o: Stock tank oil specific gravity (SG relative to water). If undefined will calculate from api
        api: Stock tank oil density (deg API). If undefined will calculate from sg_o. If both defined api value will prevail
        denomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   SWMH: Standing, Witte, McCain-Hill (1995) - Default
    """
    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    denomethod = validate_methods(["denomethod"], [denomethod])

    if sg_g == 0 and sg_sp == 0:
        print(
            "Must define at least one of sg_g and sg_sp for density calculation"
        )
        sys.exit()

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
                )  # pseudoliquid density, Eq 3.18b. Note equation in origiganl paper uses sg_sp rather than sg_g as in book.
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
        print("Must supply either sg_o or api")
        sys.exit()

    if api == 0:  # Set api from sg_o
        api = 141.5 / sg_o - 131.5
    else:  # overwrite sg_o with api value
        sg_o = oil_sg(api)

    if (
        p > pb
    ):  # Use Eq 3.20, calculating oil density from density at Pb and compressibility factor
        return fn_dic["PGTPB"](
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

    return fn_dic[denomethod.name](
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
) -> float:
    """ Returns oil formation volume factor calculated with different correlations

        p: Pressure (psia)
        pb: Bubble point pressure (psia). Defaults to 1E6, and not used for densities below Pb. A valid value is required for density calculations above Pb
        degf: Reservoir Temperature (deg F)
        rs: Oil solution gas volume (scf/stb)
        rsb: Oil solution gas volume (scf/stb) at Pb
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not known, it will be estimated from sg_sp
        sg_sp: Separator gas specific gravity (relative to air).
        sg_o: Stock tank oil specific gravity (SG relative to water).
        bomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation
                   MCAIN: McCain approach, calculating from densities
        denomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   SWMH: Standing, Witte, McCain-Hill (1995) - Default
    """
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


def oil_viso(p: float, api: float, degf: float, pb: float, rs: float) -> float:
    """ Returns Oil Viscosity with Beggs-Robinson (1975) correlation at saturated pressures
        and Petrosky-Farshad (1995) at undersaturated pressures

        p: Pressure (psia)
        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F)
        pb: Bubble point Pressure (psia)
        rs: Solution GOR (scf/stb)
    """

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


def make_bot_og(
    pi: float,
    api: float,
    degf: float,
    sg_g: float,
    pmax: float,
    pb: float = 0,
    rsb: float = 0,
    pmin: float = 25,
    nrows: int = 20,
    wt: float = 0,
    ch4_sat: float = 0,
    comethod: co_method = co_method.EXPLT,
    zmethod: z_method = z_method.DAK,
    rsmethod: rs_method = rs_method.VELAR,
    cmethod: c_method = c_method.PMC,
    denomethod: deno_method = deno_method.SWMH,
    bomethod: bo_method = bo_method.MCAIN,
    pbmethod: pb_method = pb_method.VELAR,
    export: bool = False,
    pvto: bool = False,
) -> dict:
    """
    Creates data required for Oil-Gas-Water black oil tables
    Returns dictionary of results, with index:
      - bot: Pandas table of blackoil data (for PVTO == False), or Saturated properties to pmax (if PVTO == True)
      - deno: ST Oil Density (lb/cuft)
      - deng: ST Gas Density (lb/cuft)
      - denw: Water Density at Pi (lb/cuft),
      - cw: Water Compressibility at Pi (1/psi)
      - uw: Water Viscosity at Pi (cP))
      - pb: Bubble point pressure either calculated (if only Rsb provided), or supplied by user
      - rsb: Solution GOR at Pb either calculated (if only Pb provided), or supplied by user
      - rsb_scale: The scaling factor that was needed to match user supplied Pb and Rsb
      - usat: a list of understaurated values (if PVTO == True) [usat_p, usat_bo, usat_uo]. This will be empty if PVTO == False

    If user species Pb or Rsb only, the corresponding property will be calculated
    If both Pb and Rsb are specified, then Pb calculations will be adjusted to honor both

    pi: Initial reservoir pressure (psia). Used to return water properties at initial pressure
    pb: Bubble point pressure (psia)
    rsb: Oil solution GOR at Pb (scf/stb)
    degf: Reservoir Temperature (deg F)
    sg_g: Weighted average specific gravity of surface gas (relative to air).
    api: Stock tank oil density (deg API).
    pmax: Maximum pressure to calcuate table to
    pmin: Minimum pressure to calculate table to. Default = 25
    nrows: Number of BOT rows. Default = 20
    wt: Salt wt% (0-100) in brine. Default = 0
    ch4_sat: Degree of methane saturation (0 - 1) in brine. Default = 0
    export: Boolean flag that controls whether to export full table to excel, and separate PVDG and PVDO include files. Default is False
    pvto: Boolean flag that controls whether the pvto live oil Eclipse format will be generated. This can only be active if export flag is also True;
          - extends bubble point line up to maximum pressure
          - generates undersaturated oil propeties
          - writes out PVTO include file
    """
    (
        zmethod,
        rsmethod,
        cmethod,
        denomethod,
        bomethod,
        pbmethod,
    ) = validate_methods(
        [
            "zmethod",
            "rsmethod",
            "cmethod",
            "denomethod",
            "bomethod",
            "pbmethod",
        ],
        [zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod],
    )

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=0)
    pmin = max(pmin, psc)
    sg_o = oil_sg(api)
    rsb_frac = 1.0

    # If PVTO = False
    #       - If Pb only provided, calculate rsb at Pb
    #       - If Rsb only provided, calculate Pb
    #       - If Pb AND Rsb provided, find rsb_frac (Rsb = Rsb_calc / rsb_frac) that honors both Pb and Rsb

    # If PVTO = True
    #       - If rsb_i only provided;
    #           - Calculate Pb_i at rsb_i given
    #           - Find rsb and rsb_frac at Maximum pressure point, which also delivers rs = rsb_i at Pb_i calculated
    #       - If Pb only provided;
    #           - Calculate rsb_i at Pb given
    #           - Find rsb and rsb_frac at Maximum pressure point, which also delivers rs = rsb_i at Pb given
    #       - If Pb AND Rsb provided
    #           - Find rsb and rsb_frac at Maximum pressure point, which also delivers rs = rsb_i at Pb_i calculated

    rsb_i = rsb  # _i stand for 'Initial' values given by user
    pb_i = pb
    rsb_frac = (
        1  # This is the fractional increase in rsb that Pb needs for solution.
    )

    # Calculate rsb from pb
    if rsb_i <= 0 and pb_i > 0:
        rsb = oil_rs_bub(
            degf=degf,
            api=api,
            sg_sp=sg_sp,
            sg_g = sg_g,
            pb=pb,
            rsmethod=rsmethod,
        )
        # rsb_i = rsb

    # Calculate pb from rsb
    if pb_i <= 0 and rsb_i > 0:
        pb = oil_pbub(
            degf=degf, api=api, sg_sp=sg_g, rsb=rsb, pbmethod=pbmethod
        )
        # pb_i = pb

    # Both have been defined by user. Need to work out scalar to apply to rsb to satisfy
    rsbnew = rsb
    rsb_frac = 1
    if pb_i > 0 and rsb_i > 0:
        if not pvto:  # No need to solve this if just resolving later
            #print(
            #    "Iteratively solving for Rsb fraction to use in order to harmonize user specified Pb and Rsb\n"
            #)
            pbcalc = oil_pbub(
                degf=degf, api=api, sg_sp=sg_sp, sg_g = sg_g, rsb=rsb, pbmethod=pbmethod
            )
            err = 100
            rsb_old = rsb
            i = 0
            while err > 0.0001:
                rsbnew = pb / pbcalc * rsb_old
                pbcalc = oil_pbub(
                    degf=degf,
                    api=api,
                    sg_sp=sg_sp,
                    sg_g = sg_g,
                    rsb=rsbnew,
                    pbmethod=pbmethod,
                )
                rsb_old = rsbnew
                err = abs(pb - pbcalc)
                i += 1
                if i > 100:
                    print(
                        "Could not solve Pb & Rsb for these combination of inputs"
                    )
                    sys.exit()
            rsb_frac = (
                rsb_i / rsbnew
            )  # Ratio of rsb needed to satisfy rsb defined by user at pb vs that needed to calculate Pb
    rsb_max = rsb

    if (
        pvto and pmax > pb
    ):  # PVTO has been requested, and Pb is less than max pressure requested
        # Need to find new rsb_frac that delivers rsb_i at pb_i (after depletion from pmax)
        #print(
        #    "Iteratively solving for Rsb fraction to use at maximum pressure to deliver appropriate Pb and Rsb\n"
        #)
        
        #print('degf, api, sg_sp, sg_g, pmax', degf, api, sg_sp, sg_g, pmax)
        rsb_max = oil_rs_bub(
            degf=degf,
            api=api,
            sg_sp=sg_sp,
            sg_g = sg_g,
            pb=pmax,
            rsmethod=rsmethod,
        )
        #print('rsb_max', rsb_max) #*********
        rs_at_pbi = oil_rs(
            api=api,
            degf=degf,
            sg_sp=sg_sp,
            p=pb,
            pb=pmax,
            rsb=rsb_max,
            rsmethod=rsmethod,
            pbmethod=pbmethod,
        )
        #print(rs_at_pbi)
        err = rs_at_pbi - rsb
        rsb_old = rs_at_pbi
        rsb_frac_new = rsb_frac
        i = 0
        while abs(err) > 0.0001:
            rsb_frac_new = rsb / rs_at_pbi * rsb_frac
            #print(rsb_max, rsb_frac_new, rsb_max * rsb_frac_new)
            rs_at_pbi = oil_rs(
                api=api,
                degf=degf,
                sg_sp=sg_sp,
                p=pb,
                pb=pmax,
                rsb=rsb_max * rsb_frac_new,
                rsmethod=rsmethod,
                pbmethod=pbmethod,
            )
            rsb_frac = rsb_frac_new
            err = rs_at_pbi - rsb
            i += 1
            if i > 100:
                print(
                    "Could not solve Pb & Rsb for these combination of inputs"
                )
                sys.exit()
        rsb_frac = rsb_frac_new
        
    pmax = max(pb, pmax)
    pbi = pb
    sg_sp = sg_g
    drows = 3
    if pmin in [pb, pi]:
        drows -= 1
    if pmax in [pb, pi]:
        drows -= 1
    if pb == pi:
        drows -= 1

    incr = (pmax - pmin) / (nrows - drows)

    pressures = list(pmin + incr * np.arange(0, nrows - drows + 1))
    pressures.append(pbi)
    pressures.append(pi)
    pressures = list(set(pressures))
    pressures.sort()
    pressures = np.array(pressures)
    co, cg, rss, bos, uos, gfvf, visg, gz, rvs, sg_rs, bws, visws, denos = [
        [] for x in range(13)
    ]

    if pvto:
        pb = pmax
        rsb = rsb_max * rsb_frac

    
    for p in pressures:
        if p > pb:
            rss.append(rsb)
        else:
            rss.append(
                oil_rs(
                    api=api,
                    degf=degf,
                    sg_sp=sg_sp,
                    p=p,
                    pb=pb,
                    rsb=rsb / rsb_frac,
                    rsmethod=rsmethod,
                    pbmethod=pbmethod,
                )
                * rsb_frac
            )

        bos.append(
            oil_bo(
                p=p,
                pb=pb,
                degf=degf,
                rs=rss[-1],
                rsb=rsb,
                sg_g=sg_g,
                sg_sp=sg_sp,
                sg_o=sg_o,
                denomethod=denomethod,
                bomethod=bomethod,
            )
        )
        denos.append(
            oil_deno(
                p=p,
                degf=degf,
                rs=rss[-1],
                rsb=rsb,
                sg_g=sg_g,
                sg_sp=sg_sp,
                pb=pb,
                sg_o=sg_o,
                api=api,
            )
        )
        uos.append(oil_viso(p=p, api=api, degf=degf, pb=pb, rs=rss[-1]))
        co.append(
            oil_co(
                p=p,
                api=api,
                sg_sp=sg_sp,
                sg_g=sg_g,
                degf=degf,
                pb=pb,
                rsb=rss[-1],
                zmethod=zmethod,
                rsmethod=rsmethod,
                cmethod=cmethod,
                denomethod=denomethod,
                bomethod=bomethod,
            )
        )

        gfvf.append(
            gas.gas_bg(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
            * 1000
            / CUFTperBBL
        )  # rb/mscf
        gz.append(
            gas.gas_z(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        visg.append(
            gas.gas_ug(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        cg.append(gas.gas_cg(p=p, sg=sg_g, degf=degf, cmethod=cmethod))
        bw, lden, visw, cw, rsw = brine.brine_props(
            p=p, degf=degf, wt=wt, ch4_sat=ch4_sat
        )
        bws.append(bw)
        visws.append(visw)

    # And undersaturated lines if required
    if pvto:
        usat_bo = []
        usat_uo = []
        usat_p = []

        for i, p in enumerate(pressures):
            if i == 0:
                continue
            try:
                usat_p.append(pressures[i:])
                usat_bo.append(
                    [
                        oil_bo(
                            p=pusat,
                            pb=p,
                            degf=degf,
                            rs=rss[i],
                            rsb=rss[i],
                            sg_g=sg_g,
                            sg_sp=sg_sp,
                            sg_o=sg_o,
                            denomethod=denomethod,
                            bomethod=bomethod,
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
            except:
                pass

    st_deno = sg_o * WDEN  # lb/cuft
    st_deng = gas.gas_den(
        p=psc, sg=sg_g, degf=tsc, zmethod=zmethod, cmethod=cmethod
    )
    bw, lden, visw, cw, rsw = brine.brine_props(
        p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat
    )
    res_denw = lden * WDEN  # lb/cuft
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

        if pvto:  # Also export PVTO include file;
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
                except:
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
