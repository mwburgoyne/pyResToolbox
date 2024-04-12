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
from scipy.integrate import quad
from scipy.optimize import brentq
#from scipy.optimize import minimize
import pandas as pd
from tabulate import tabulate
from gwr_inversion import gwr
from mpmath import mp

#import simtools.simtools as simtools
from .simtools import simtools 

# Constants
R = 10.731577089016  # Universal gas constant, ft³·psia/°R·lb.mol
psc = 14.696  # Standard conditions pressure (psia)
tscf = 60  # Standard conditions temperature (deg F)
f2r = 459.67  # Offset to convert degrees F to degrees Rankine
tscr = tscf + f2r  # Standard conditions temperature (deg R)
mw_air = 28.97  # MW of Air
scf_per_mol = R * tscr / psc  # scf/lb-mol (V = ZnRT/P, Z = 1, n = 1)


class z_method(Enum):  # Gas Z-Factor calculation model
    DAK = 0
    HY = 1
    WYW = 2
    PR = 3


class c_method(Enum):  # Gas critical properties calculation method
    PMC = 0
    SUT = 1


class pb_method(Enum):  # Bubble point calculation method
    STAN = 0
    VALMC = 1
    VELAR = 2


class rs_method(Enum):  # Oil solution gas calculation method
    VELAR = 0
    STAN = 1
    VASBG = 2


class bo_method(Enum):  # Oil FVF calculation method
    MCAIN = 0
    STAN = 1


class uo_method(Enum):  # Oil viscosity calculation method
    BR = 0


class deno_method(Enum):  # Oil Density calculation method
    SWMH = 0


class co_method(Enum):  # Oil compressibility calculation method
    EXPLT = 0


class kr_family(Enum):  # Relative permeability family type
    COR = 0
    LET = 1


class kr_table(Enum):  # Relative permeability table type
    SWOF = 0
    SGOF = 1
    SGWFN = 2


class_dic = {
    "zmethod": z_method,
    "cmethod": c_method,
    "pbmethod": pb_method,
    "rsmethod": rs_method,
    "bomethod": bo_method,
    "uomethod": uo_method,
    "denomethod": deno_method,
    "comethod": co_method,
    "krfamily": kr_family,
    "krtable": kr_table,
}


def bisect_solve(args, f, xmin, xmax, rtol):
    err_hi = f(args, xmax)
    err_lo = f(args, xmin)
    iternum = 0
    err_mid = 1
    while abs(err_mid) > rtol:
        mid_val = (xmax + xmin) / 2
        err_mid = f(args, mid_val)
        iternum += 1
        if iternum > 99:
            print("Could not solve via bisection")
            sys.exit()
        if (err_hi * err_mid < 0):  # Solution point must be higher than current mid_val case
            xmin = mid_val
            err_lo = err_mid
            mid_val = (mid_val + xmax) / 2
        else:
            xmax = mid_val  # Other_wise must be lower than current mid_val case
            err_hi = err_mid
    return mid_val


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
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
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
                 'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified
        cmethod: Method for calculating critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               Defaults to 'PMC'
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        S: Skin. Defaults to zero if undefined
        D: Non Darcy Skin Factor (day/mscf). Defaults to zero if undefined
        sg: Gas SG relative to air, Defaults to 0.75 if undefined

    """
    k, h, pr, pwf = (
        np.asarray(k),
        np.asarray(h),
        np.asarray(pr),
        np.asarray(pwf),
    )
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    tc, pc = gas_tc_pc(sg, n2, co2, h2s, cmethod.name, tc, pc)

    direction = 1
    if pr.size + pwf.size == 2:  # Single set of pressures
        if pr < pwf:
            direction = -1  # Direction is needed because solving the quadratic with non-Darcy factor will fail if using a negative delta_mp
        delta_mp = abs(
            gas_dmp(
                p1=pwf,
                p2=pr,
                degf=degf,
                sg=sg,
                zmethod=zmethod,
                cmethod=cmethod,
                tc=tc,
                pc=pc,
                n2=n2,
                co2=co2,
                h2s=h2s,
            )
        )
    else:
        if pr.size > 1:  # Multiple Pr's
            direction = np.array([p > pwf for p in pr])
            direction = 2 * direction - 1
            delta_mp = np.absolute(
                np.array(
                    [
                        gas_dmp(
                            p1=p,
                            p2=pwf,
                            degf=degf,
                            sg=sg,
                            zmethod=zmethod,
                            cmethod=cmethod,
                            tc=tc,
                            pc=pc,
                            n2=n2,
                            co2=co2,
                            h2s=h2s,
                        )
                        for p in pr
                    ]
                )
            )
        else:  # Multiple BHFP's
            direction = np.array([pr > bhfp for bhfp in pwf])
            direction = 2 * direction - 1
            delta_mp = np.absolute(
                np.array(
                    [
                        gas_dmp(
                            p1=pr,
                            p2=bhfp,
                            degf=degf,
                            sg=sg,
                            zmethod=zmethod,
                            cmethod=cmethod,
                            tc=tc,
                            pc=pc,
                            n2=n2,
                            co2=co2,
                            h2s=h2s,
                        )
                        for bhfp in pwf
                    ]
                )
            )

    qg = darcy_gas(delta_mp, k, h, degf, r_w, r_ext, S, D, radial=True)
    return direction * qg


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
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
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
                 'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified
        cmethod: Method for calculting critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               Defaults to 'PMC' if not specified
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if not specified
        co2: Molar fraction of CO2. Defaults to zero if not specified
        h2s: Molar fraction of H2S. Defaults to zero if not specified
        sg: Gas SG relative to air, Defaults to 0.75 if not specified
        degf: Reservoir Temperature (deg F).
    """
    k, area, pr, pwf = (
        np.asarray(k),
        np.asarray(area),
        np.asarray(pr),
        np.asarray(pwf),
    )
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    tc, pc = gas_tc_pc(sg, n2, co2, h2s, cmethod.name, tc, pc)

    direction = 1
    if pr.size + pwf.size == 2:  # Single set of pressures
        if pr < pwf:
            direction = (
                -1
            )  # Direction is needed because solving the quadratic with non-Darcy factor will fail if using a negative delta_mp
        delta_mp = abs(
            gas_dmp(
                p1=pwf,
                p2=pr,
                degf=degf,
                sg=sg,
                zmethod=zmethod,
                cmethod=cmethod,
                tc=tc,
                pc=pc,
                n2=n2,
                co2=co2,
                h2s=h2s,
            )
        )
    else:
        if pr.size > 1:
            direction = np.array([p > pwf for p in pr])
            direction = 2 * direction - 1
            delta_mp = np.absolute(
                np.array(
                    [
                        gas_dmp(
                            p1=p,
                            p2=pwf,
                            degf=degf,
                            sg=sg,
                            zmethod=zmethod,
                            cmethod=cmethod,
                            tc=tc,
                            pc=pc,
                            n2=n2,
                            co2=co2,
                            h2s=h2s,
                        )
                        for p in pr
                    ]
                )
            )
        else:
            direction = np.array([pr > bhfp for bhfp in pwf])
            direction = 2 * direction - 1
            delta_mp = np.absolute(
                np.array(
                    [
                        gas_dmp(
                            p1=pr,
                            p2=bhfp,
                            degf=degf,
                            sg=sg,
                            zmethod=zmethod,
                            cmethod=cmethod,
                            tc=tc,
                            pc=pc,
                            n2=n2,
                            co2=co2,
                            h2s=h2s,
                        )
                        for bhfp in pwf
                    ]
                )
            )

    qg = darcy_gas(delta_mp, k, 1, degf, area, length, 0, 0, radial=False)
    return direction * qg


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
    tr = degf + f2r
    if radial:
        a = k * h * delta_mp
        b = 1422 * tr
        c = np.log(l2 / l1) - 0.75 + S
        if D > 1e-9:  # Solve analytically for rate with non-Darcy factor by rearranging into root of a quadratic equation.
            return (np.sqrt(4 * a * b * D + (b * b * c * c)) - (b * c)) / (2 * b * D)
    else:
        a = k * h * l1 * delta_mp
        b = 2 * np.pi * 1422 * tr
        c = l2
    # Else, ignore non-Darcy skin
    return a / (b * c)


def gas_tc_pc(
    sg: float,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    cmethod: str = "PMC",
    tc: float = 0,
    pc: float = 0,
) -> tuple:
    """ Returns a tuple of critical temperature (deg R) and critical pressure (psia) for hydrocarbon gas
        cmethod: 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        sg: Specific gravity of reservoir gas (relative to air)
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
    """
    if tc * pc > 0:  # Critical properties have both been user specified
        return (tc, pc)
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

    elif (
        cmethod.name == "SUT"
    ):  # Sutton equations with Wichert & Aziz corrections
        sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / mw_air) / (
            1 - n2 - co2 - h2s
        )  # Eq 3.53
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

        # Changed to SBV mixing rules per SPE 14265
        # yis = np.array([1-n2-co2-h2s, n2, co2, h2s])
        # if yis[0] == 1.0:
        #    ppc = ppc_hc
        #    tpc = tpc_hc
        # else:
        #    tcs = np.array([tpc_hc, 239.26, 547.58, 672.35])
        #    pcs = np.array([ppc_hc, 507.5, 1071.0, 1306.0])
        #    J = np.sum(yis*tcs/pcs)/3 + 2*(np.sum(yis*np.sqrt(tcs/pcs)))**2/3  # Eq 16
        #    K = np.sum(yis*tcs/pcs**0.5)         # Eq 17
        #    tpc = K**2/J                         # Eq 18
        #    ppc = tpc/J                          # Eq 19

    else:
        print("Incorrect cmethod specified")
        sys.exit()

    if tc > 0:
        tpc = tc
    if pc > 0:
        ppc = pc
    return (tpc, ppc)


def gas_z(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
) -> np.ndarray:
    """ Returns real-gas deviation factor (Z). Returning either single float, or numpy array depending upon 
        whether single pressure of list/array or pressures has been specified.
        p: Gas pressure (psia). Takes a single float, 1D list or 1D Numpy array
        sg: Gas SG relative to air. Single float only
        degf: Gas Temperature (deg F). Single float only
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified
        cmethod: Method for calculting critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               Defaults to 'PMC'
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
    """
    p = np.asarray(p)
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    tc, pc = gas_tc_pc(sg, n2, co2, h2s, cmethod.name, tc, pc)
    tr = (degf + f2r) / tc
    
    single_p = False
    if p.size == 1:
        single_p = True
        pprs = np.array([p/pc])
    else:
        pprs = np.array(p/pc)
    
    def zdak(pprs, tr):
        # DAK from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        # sg relative to air, t in deg F, p in psia, n2, co2 and h2s in fractions (0-1)
        
        tr2, tr3 = tr**2, tr**3
        a = np.array([0,0.3265,-1.07,-0.5339,0.01569,-0.05165,0.5475,-0.7361,0.1844,0.1056,0.6134,0.7210])
        c1 = a[1]+(a[2]/tr) + (a[3]/tr3) + (a[4]/tr**4) + (a[5]/tr**5)
        c2 = a[6] + (a[7]/tr) + (a[8]/tr2)
        c3 = a[9] * ((a[7]/tr) + (a[8]/tr2))
    
        # Unpublished empirical fit of reduced pressure at minimum DAK-Z as a function of 
        # reduced Temperature over Tpr range 1.03 - 4.0
        # By Mark Burgoyne, July 2023
        def min_ppr(tpr): #  N. 127:   Y = A/X+B*EXP(C*X)+D <-- Hyp + exp + c
            if tpr > 3.43:
                return 0.0
            if tpr < 1.03:
                return 1.2625075898234461
            A, B, C, D = 0.2283E+02, -.1359E+03, -.2188E+01, -.6614E+01
            minppr = A/tpr+B*np.exp(C*tpr)+D
            return minppr
            
        # Unpublished empirical fit of minimum DAK-Z values as a function of 
        # reduced Temperature over Tpr range 1.03 - 4.0
        # By Mark Burgoyne, July 2023
        def fit_minz(tpr): # N. 149:   Y = (A+X)/(B+C*X**2)+D <--- Straight line / parabola + c
            if tpr > 3.43:
                return 1.0
            A, B, C, D = -.1556E+01, 0.8930E-01, 0.8043E+00, 0.8019E+00
            minzs = (A+tpr)/(B+C*tpr**2)+D 
            return minzs
    
        def new_dak_z(z, pr, tr):
            rhor = 0.27 * pr / (tr * z)  # 2.8
            c4 = a[10] * (1 + a[11] * rhor**2) * (rhor**2/tr3) * np.exp(-a[11] * rhor**2)
            err = fz(z, rhor, c4)
            return z - err/fz_prime(z, rhor, c4), err # Returns updated guess for Z, as well as the error
        
        def fz(z, rhor, c4):                          # The DAK Error function
            return z - (1 + c1*rhor + c2*rhor**2 - c3 * rhor**5 + c4)
    
        def fz_prime(z, rhor, c4):                    # Derivative of the DAK Error function
            return 1 + c1*rhor/z + (2*c2*rhor**2/z) - (5*c3*rhor**5/z) + 2*a[10]*rhor**2/(z*tr3)*(1+a[11]*rhor**2 - (a[11]*rhor**2)**2)*np.exp(-a[11]*rhor**2)
        
        def z_err(z, *args):
            pr = args[0]
            rhor = 0.27 * pr / (tr * z)  # 2.8
            c4 = a[10] * (1 + a[11] * rhor**2) * (rhor**2/tr3) * np.exp(-a[11] * rhor**2)
            z2 = (1 + c1*rhor + c2*rhor**2 - c3 * rhor**5 + c4)
            return z2 - z        
            
        def z_dak_calc(pr, tr, tol = 1e-6):
            niter = 0
            minppr = min_ppr(tr)
            
            if abs(pr - minppr) > 0.05: # If Ppr is further from calculated minimum Ppr than 0.05, use Newton solver
                newton_solve = True
            else:
                newton_solve = False    # Else, use bisection solver, and setup Z bounds to search within
                minz = fit_minz(tr)
                bounds = (minz - 0.02, minz + 0.02)
    
            if newton_solve:
                midz = min(max(0.1, z_pr([pr], tr)),3) # First guess using explicit calculation method
                mid_err = 1
                while abs(mid_err) > tol and niter < 100:
                    midz, mid_err = new_dak_z(midz, pr, tr)
                    niter += 1
            else:
                midz = brentq(z_err, bounds[0], bounds[1], args=(pr))  
            return midz
        
        zout = [z_dak_calc(pr, tr) for pr in pprs]
            
        if single_p:
            return zout[0]
        else:
            return np.array(zout)

    # Hall & Yarborough
    def z_hy(pprs, tr):
        tr = 1/tr
        t2 = tr ** 2
        a = 0.06125 * tr * np.exp(-1.2 * (1 - tr) ** 2)
        b = tr * (14.76 - 9.76 * tr + 4.58 * t2)
        c = tr * (90.7 - 242.2 * tr + 42.4 * t2)
        D = 2.18 + 2.82 * tr
         
        # f(y)
        def f(y, a, b, c, D, pr):   
            return ((y + y ** 2 + y ** 3 - y ** 4) / ((1 - y) ** 3)) - a * pr - b * y ** 2 + c * y ** D
    
        # derivative of f(y)
        def df(y, a, b, c, D):
            return ((1 + 4 * y + 4 * y ** 2 - 4 * y ** 3 + y ** 4) / ((1 - y) ** 4)) - 2 * b * y + c * D * y ** (D - 1)     
         
        zout = []
        for pr in pprs:
            yi = a*pr / z_wyw(pr, 1/tr) # First guess
            niter, y = 0, 0.01
            while (abs(y-yi)/y) > 0.0005 and niter < 100: 
                # Newton Raphson
                y = yi - (f(yi, a, b, c, D, pr) / df(yi, a, b, c, D))
                niter += 1
                yi = y
            zout.append(a * pr / y)
                
        if single_p:
            return zout[0]
        else:
            return np.array(zout)

    # Wang, Ye & Wu, 2021, 0.2 < Ppr < 30, 1.05 < tpr < 3.0
    # "An accurate correlation for calculating natural gas compressibility factors under a wide range of pressure conditions"
    # https://doi.org/10.1016/j.egyr.2021.11.029
    def z_wyw(pprs, tr):
        a = [0, 256.41675, 7.18202, -178.5725, 182.98704, -40.74427, 2.24427, 47.44825, 5.2852, -0.14914, 271.50446, 16.2694, -121.51728, 167.71477, -81.73093, 20.36191, -2.1177, 124.64444, -6.74331, 0.20897, -0.00314]
        numerators = a[1] + a[2] * (1 + a[3] * tr + a[4] * tr ** 2 + a[5] * tr ** 3 + a[6] * tr ** 4) * pprs + a[7] * pprs ** 2 + a[8] * pprs ** 3 + a[9] * pprs ** 4
        denominators = a[10] + a[11] * (1 + a[12] * tr + a[13] * tr ** 2 + a[14] * tr ** 3 + a[15] * tr ** 4 + a[16] * tr ** 5) * pprs + a[17] * pprs ** 2 + a[18] * pprs ** 3 + a[19] * pprs ** 4 + a[20] * pprs ** 5
        zs = numerators/denominators

        if single_p:
            return float(zs)
        else:
            return zs

    def z_pr(pprs, tr):
        
        # Analytic solution for maximum real root of cubic polynomial
        # a[0] * Z^3 + a[1]*Z^2 + a[2]*Z + a[3] = 0
        def max_root(a):
            p = (3 * a[2]- a[1]**2)/3
            q = (2 * a[1]**3 - 9 * a[1] * a[2] + 27 * a[3])/27
            root_diagnostic = q**2/4 + p**3/27

            if root_diagnostic < 0:
                m = 2*np.sqrt(-p/3)
                qpm = 3*q/p/m
                theta1 = np.arccos(qpm)/3
                roots = np.array([m*np.cos(theta1), m*np.cos(theta1+4*np.pi/3), m*np.cos(theta1+2*np.pi/3)])
                Zs = roots - a[1] / 3
            else:
                P = (-q/2 + np.sqrt(root_diagnostic))
                if P >= 0:
                    P = P **(1/3)
                else:
                    P = -(-P)**(1/3)
            
                Q = (-q/2 - np.sqrt(root_diagnostic))
                if Q >=0:
                    Q = Q **(1/3)
                else:
                    Q = -(-Q)**(1/3)
                Zs = np.array([P + Q]) - a[1] / 3
            return max(Zs)
    
        w, s = -0.048964, -0.394899
        m = 0.37464 + 1.54226 * w - 0.26992 * w**2
        alpha = (1 + m * (1 - np.sqrt(tr)))**2    
        
        zout = []
        for pr in pprs:
            A = 0.429188 * alpha * pr / tr**2
            B = 0.0692551 * pr / tr
    
            # Coefficients of Cubic: a[0] * Z^3 + a[1]*Z^2 + a[2]*Z + a[3] = 0
            a = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]  
            zout.append(max_root(a) - s * B) # Volume translated Z          
                
        if single_p:
            return zout[0]
        else:
            return np.array(zout)
            
        if single_p:
            return float(zs)
        else:
            return zs
            
    zfuncs = {"DAK": zdak, "HY": z_hy, "WYW": z_wyw, "PR": z_pr}

    return zfuncs[zmethod.name](pprs, tr)


def gas_ug(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
) -> np.ndarray:
    """ Returns Gas Viscosity (cP)
        Uses Lee, Gonzalez & Eakin (1966) Correlation using equations 2.14-2.17 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.

          p: Gas pressure (psia)
          sg: Gas SG relative to air
          degf: Reservoir Temperature (deg F).
          zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                   defaults to 'DAK' if not specified
          cmethod: Method for calculting critical properties
                   'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                   'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   Defaults to 'PMC'
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
    """
    p = np.asarray(p)
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    zee = gas_z(
        p=p,
        sg=sg,
        degf=degf,
        zmethod=zmethod,
        cmethod=cmethod,
        tc=tc,
        pc=pc,
        n2=n2,
        co2=co2,
        h2s=h2s,
    )
    t = degf + f2r
    m = mw_air * sg
    rho = m * p / (t * zee * R * 62.37)
    b = 3.448 + (986.4 / t) + (0.01009 * m)  # 2.16
    c = 2.447 - (0.2224 * b)  # 2.17
    a = (
        (9.379 + (0.01607 * m)) * np.power(t, 1.5) / (209.2 + (19.26 * m) + t)
    )  # 2.15
    return a * 0.0001 * np.exp(b * np.power(rho, c))  # 2.14


def gas_ugz(
    p: npt.ArrayLike, sg: float, degf: float, zee: npt.ArrayLike
) -> np.ndarray:
    """ Returns product of Gas Viscosity (cP) * Gas Z-Factor
        Uses Lee, Gonzalez & Eakin (1966) Correlation using equations 2.14-2.17 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        Same as ug function, but with a precalculated Z factor to eliminate duplicated evaluation in m(p) calculations
        p: Gas pressure (psia)
        degf: Gas Temperature (deg F)
        sg: Specific gravity of reservoir gas (relative to air)
        zee: pre-calculated gas Z-Factor
    """
    p, zee = np.asarray(p), np.asarray(zee)
    try:
        if len(p) != len(zee):
            print(
                "Warning, length of pressure and z-factor arrays should be the same"
            )
    except:
        pass
    t = degf + f2r
    m = mw_air * sg
    rho = m * p / (t * zee * R * 62.37)
    b = 3.448 + (986.4 / t) + (0.01009 * m)  # 2.16
    c = 2.447 - (0.2224 * b)  # 2.17
    a = (
        (9.379 + (0.01607 * m)) * np.power(t, 1.5) / (209.2 + (19.26 * m) + t)
    )  # 2.15
    return a * 0.0001 * np.exp(b * np.power(rho, c)) * zee  # 2.14


def gas_cg(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
    cmethod: c_method = c_method.PMC,
) -> np.ndarray:
    """ Returns gas compressibility (1/psi) using the 'DAK' Dranchuk & Abou-Kassem (1975) Z-Factor &
        Critical property correlation values if not explicitly specified
        p: Gas pressure (psia)
        sg: Gas SG relative to air. Defaults to False if undefined
        pwf: BHFP (psia)
        degf: Gas Temperature (deg F)
        cmethod: Method for calculating Tc and Pc

        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
    """
    cmethod = validate_methods(["cmethod"], [cmethod])
    p = np.asarray(p)
    tc, pc = gas_tc_pc(
        sg=sg, n2=n2, co2=co2, h2s=h2s, tc=tc, pc=pc, cmethod=cmethod
    )
    pr = p / pc
    tr = (degf + f2r) / tc
    zee = gas_z(p=p, degf=degf, sg=sg, tc=tc, pc=pc, n2=n2, co2=co2, h2s=h2s)

    a = [0, 0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210]
    rhor = 0.27 * pr / (tr * zee)
    dzdrho = (
        a[1]
        + (a[2] / tr)
        + (a[3] / (tr * tr * tr))
        + (a[4] / (tr * tr * tr * tr))
        + (a[5] / np.power(tr, 5))
    )
    dzdrho = dzdrho + (2 * rhor * (a[6] + (a[7] / tr) + (a[8] / (tr * tr))))
    dzdrho = dzdrho - (
        5 * np.power(rhor, 4) * a[9] * ((a[7] / tr) + (a[8] / (tr * tr)))
    )
    dzdrho = dzdrho + (2 * a[10] * rhor / (tr * tr * tr)) * (
        1 + (a[11] * rhor * rhor) - (a[11] * a[11] * np.power(rhor, 4))
    ) * np.exp(
        -a[11] * rhor * rhor
    )  # 2.23
    cpr = (1 / pr) - (
        (0.27 / (zee * zee * tr)) * (dzdrho / (1 + (rhor / zee) * dzdrho))
    )  # 2.22
    cg = cpr / pc  # 2.21
    return cg


def gas_bg(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
) -> np.ndarray:
    """ Returns Bg (gas formation volume factor) for natural gas (rcf/scf)
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 Defaults to 'PMC'
          p: Gas pressure (psia)
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          sg: Gas SG relative to air
          degf: Reservoir Temperature (deg F)
    """
    p = np.asarray(p)
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    zee = gas_z(
        p=p,
        degf=degf,
        sg=sg,
        tc=tc,
        pc=pc,
        n2=n2,
        co2=co2,
        h2s=h2s,
        zmethod=zmethod,
        cmethod=cmethod,
    )
    return zee * (degf + f2r) / (p * 35.37)


def gas_den(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
) -> np.ndarray:
    """ Returns gas density for natural gas (lb/cuft)

          zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                   defaults to 'DAK' if not specified
          cmethod: Method for calculting critical properties
                   'sut' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                   'pmc' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   Defaults to 'pmc'
          p: Gas pressure (psia)
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          sg: Gas SG relative to air, Defaults to False if undefined
          degf: Reservoir Temperature (deg F). Defaults to False if undefined
    """
    p = np.asarray(p)
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    zee = gas_z(
        p=p,
        degf=degf,
        sg=sg,
        tc=tc,
        pc=pc,
        n2=n2,
        co2=co2,
        h2s=h2s,
        zmethod=zmethod,
        cmethod=cmethod,
    )
    m = sg * mw_air
    t = degf + f2r
    rhog = p * m / (zee * R * t)
    return rhog


def gas_ponz2p(
    poverz: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
    rtol: float = 1e-7,
) -> np.ndarray:
    """ Returns pressure corresponding to a P/Z value for natural gas (psia)
        Calculated through iterative solution method
        poverz: Gas pressure / Z-Factor (psia)

        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 Defaults to 'PMC'
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          sg: Gas SG relative to air, Defaults to False if undefined
          degf: Reservoir Temperature (deg F). Defaults to False if undefined
          rtol: Relative solution tolerance. Will iterate until abs[(poverz - calculation)/poverz] < rtol
    """

    def PonZ2P_err(args, p):
        ponz, sg, t, zmethod, cmethod, tc, pc, n2, co2, h2s = args
        zee = gas_z(
            p=p,
            degf=t,
            sg=sg,
            tc=tc,
            pc=pc,
            n2=n2,
            co2=co2,
            h2s=h2s,
            zmethod=zmethod,
            cmethod=cmethod,
        )
        return (p - (ponz * zee)) / p

    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    poverz = np.asarray(poverz)
    single_p = False
    if poverz.size == 1:
        single_p = True
        poverz = [poverz]
    else:
        poverz = poverz.tolist()

    p = []
    for ponz in poverz:
        args = (ponz, sg, degf, zmethod, cmethod, tc, pc, n2, co2, h2s)
        p.append(bisect_solve(args, PonZ2P_err, ponz * 0.2, ponz * 1.8, rtol))
    p = np.array(p)
    if single_p:
        p = p[0]
    return p


def gas_grad2sg(
    grad: float,
    p: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
    rtol: float = 1e-7,
) -> float:
    """ Returns insitu gas specific gravity consistent with observed gas gradient. Solution iteratively calculated via bisection
        Calculated through iterative solution method. Will fail if gas SG is below 0.55, or greater than 1.75.

        grad: Observed gas gradient (psi/ft)
        p: Pressure at observation (psia)
        degf: Reservoir Temperature (deg F). Defaults to False if undefined
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 Defaults to 'PMC'
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          sg: Gas SG relative to air, Defaults to False if undefined
          rtol: Relative solution tolerance. Will iterate until abs[(grad - calculation)/grad] < rtol
    """

    t = degf + f2r

    def grad_err(args, sg):
        grad, p, zmethod, cmethod, tc, pc, n2, co2, h2s = args
        m = sg * mw_air
        zee = gas_z(
            p=p,
            degf=degf,
            sg=sg,
            tc=tc,
            pc=pc,
            n2=n2,
            co2=co2,
            h2s=h2s,
            zmethod=zmethod,
            cmethod=cmethod,
        )
        grad_calc = p * m / (zee * R * t) / 144
        error = (grad - grad_calc) / grad
        return error

    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    args = (grad, p, zmethod, cmethod, tc, pc, n2, co2, h2s)
    return bisect_solve(args, grad_err, 0.55, 1.75, rtol)


def gas_dmp(
    p1: float,
    p2: float,
    degf: float,
    sg: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    n2: float = 0,
    co2: float = 0,
    h2s: float = 0,
    tc: float = 0,
    pc: float = 0,
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
                   'PR' Tuned single component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                   defaults to 'DAK' if not specified
        cmethod: Method for calculting critical properties
                 'sut' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'pmc' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 Defaults to 'pmc'
        tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
        pc: Critical gas pressure (psia). Calculates using cmethod if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
    """

    def m_p(p, *args):
        # Pseudo pressure function to be integrated
        degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s = args
        zee = gas_z(
            p=p,
            degf=degf,
            sg=sg,
            zmethod=zmethod,
            cmethod=cmethod,
            n2=n2,
            co2=co2,
            h2s=h2s,
            tc=tc,
            pc=pc,
        )
        mugz = gas_ugz(
            p, degf, sg, zee
        )  # Gas viscosity z-factor product using a precalculated Z factor
        return 2 * p / (mugz)

    if p1 == p2:
        return 0
    zmethod, cmethod = validate_methods(
        ["zmethod", "cmethod"], [zmethod, cmethod]
    )

    return quad(
        m_p,
        p1,
        p2,
        args=(degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s),
        limit=500,
    )[0]


def gas_fws_sg(sg_g: float, cgr: float, api_st: float) -> float:
    """
     Estimates FWS specific gravity of gas-condensate from separator gas SG, CGR and API
     Uses Standing correlation to estimate condensate MW from API.
     Returns SG of FWS gas

     sg_g: Specific gravity of weighted average surface gas (relative to air)
     api_st: Density of stock tank liquid (API)
     cgr: Condensate gas ratio (stb/mmscf)
     """
    # 1 mmscf separator gas basis with 379.482 scf/lb-mole
    cond_vol = cgr * 5.61458  # cuft/mmscf surface gas
    cond_sg = oil_sg(api_st)
    cond_mass = cond_vol * (cond_sg * 62.4)  # lb
    surface_gas_moles = 1e6 / scf_per_mol  # lb-moles
    surface_gas_mass = sg_g * mw_air * surface_gas_moles  # lb
    cond_mw = 240 - 2.22 * api_st  # lb/lb-moles (Standing correlation)
    cond_moles = cond_mass / cond_mw  # lb-moles
    fws_gas_mass = cond_mass + surface_gas_mass  # lb
    fws_gas_moles = cond_moles + surface_gas_moles  # lb-moles
    fws_gas_mw = fws_gas_mass / fws_gas_moles  # lb/lb-moles
    return fws_gas_mw / mw_air


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
        S: Wellbore Skin (Dimensionless). Defaults to zero if not specified
        uo: Liquid viscosity (cP)
        bo: Liquid Formation Volume Factor (rb/stb)
        pb: Bubble point pressure (psia). Defaults to zero if not specified. Not used unless Vogel option is invoked
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
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
        pb: Bubble point pressure (psia). Defaults to zero if not specified. Not used unless Vogel option is invoked
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
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
        mw: Molecular weight of the liquid (g/gmole / lb/lb.mol)
        Ja: Varies between 0 (Paraffins) - 1 (Aromatic)
    """
    ja = min(1, ja)
    ja = max(0, ja)
    return 0.8468 - 15.8 / mw + ja * (0.2456 - 1.77 / mw)


def oil_twu_props(
    mw: float, ja: float = 0, sg: float = 0, damp: float = 1
) -> tuple:
    """ Returns tuple of tb, tc, vc, pc using method from Twu (1984) correlations for petroleum liquids
        Modified with damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu)
        Returns sg, tb (R), tc (R), pc (psia), vc (ft3/lbmol)

        mw: Molecular weight of the liquid hydrocarbon (g/g.mol / lb/lb.mol)
        ja: jacoby Aromaticity Factor relationship. Varies between 0 (Paraffins) - 1 (Aromatic). Defaults to zero if undefined
        sg: Specific gravity of the liquid (fraction relative to water density). Will use jacoby method to estimate sg from mw if undefined.
        damp: damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu). Defaults to 1
        Unless otherwise mentioned, all Twu equation references
        are from Whitson Monograph
    """
    if sg == 0:
        sg = oil_ja_sg(
            mw, ja
        )  # Use jacoby relationship to estimate sg if not specified

    # Estimate boiling point given mw, sg and Paraffinicity
    # damp = 0 (Paraffin) - 1 (Original Twu)
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
                print("Check inputs. Twu algorithm did not converge")
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
    tc = Twu_tc(tb, sgp, sg)
    vc = Twu_vc(tb, tcp, sg, sgp)
    pc = Twu_pc(tb, sgp, sg, pcp, tc, tcp, vc, vcp)
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
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
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
            print(rsb, api, sg_sp, degf)
            print(
                "Need valid values for rsb, api, sg_sp and degf for Velarde or Valko McCain Pb calculation"
            )
            sys.exit()

    def pbub_standing(
        api, degf, sg_g, rsb, sg_sp
    ) -> float:  # 1.63 in 'Oil 7 Gas Properties & Correlations' - http://dx.doi.org/10.1016/B978-0-12-803437-8.00001-4
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
    pbmethod: pb_method = pb_method.VALMC,
    rsmethod: rs_method = rs_method.VELAR,
) -> float:
    """ Returns Solution GOR (scf/stb) at bubble point pressure.
        Uses the inverse of the Bubble point pressure correlations, with the same method families
        Note: At low pressures, the VALMC method will fail (generally when Rsb < 10 scf/stb).
              The VALMC method will revert to the STAN method in these cases

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F)
        pb: Bubble point Pressure (psia)
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - Default
                   VELAR: Velarde, Blasingame & McCain (1997)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
    """
    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )

    if pbmethod.name == "STAN":
        if pb * api * sg_g * degf == 0:
            print(
                "Need valid values for pb, api, sg_g for Standing Correlation"
            )
            sys.exit()
    else:
        if pb * api * sg_sp * degf == 0:
            print(
                "Need valid values for pb, api, sg_sp and degf for Velarde or Valko McCain Pb calculation"
            )
            sys.exit()

    def rsbub_standing(api, degf, pb, sg_g, sg_sp) -> float:
        a = 0.00091 * degf - 0.0125 * api  # Eq 1.64
        return sg_sp * (((pb - psc) / 18.2 + 1.4) / 10 ** a) ** (
            1 / 0.83
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def rsbub_valko_mccain(api, degf, pb, sg_g, sg_sp) -> float:
        # Solve via iteration. First guess using Standing Rsb, then simple Newton Iterations
        old_rsb = rsbub_standing(api, degf, pb, sg_g, sg_sp)
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
            error_slope = (new_err - old_err) / (new_rsb - old_rsb)
            old_err, old_rsb = new_err, new_rsb

            intcpt = new_err - error_slope * new_rsb

            new_rsb = -intcpt / error_slope
            if (
                i > 100
            ):  # At low rsb VALMC will not converge, use Velarde instead
                return rsbub_velarde(api, degf, pb, sg_g, sg_sp)
        return new_rsb

    def rsbub_velarde(api, degf, pb, sg_g, sg_sp) -> float:
        x = 0.013098 * degf ** 0.282372 - 8.2e-6 * api ** 2.176124
        rsb = (
            0.270811 * sg_sp ** (10093 / 62500) * pb ** 0.186745 * 10 ** (-x)
            + 92519 * sg_sp ** (10093 / 62500) * 2 ** (-x - 3) * 5 ** (-x - 6)
        ) ** (200000 / 16293)
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


def validate_methods(names, variables):
    for m, method in enumerate(names):
        if type(variables[m]) == str:
            try:
                variables[m] = class_dic[method][variables[m].upper()]
            except:
                print("An incorrect method was specified")
                sys.exit()
    if len(variables) == 1:
        return variables[0]
    else:
        return variables


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

    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
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
            pbmethod=pbmethod,
            rsmethod=rsmethod,
        )

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
                * np.exp(25.7240 * (api / (degf + f2r)))
            )
        else:
            return (
                0.0178
                * sg_gs
                * p ** 1.1870
                * np.exp(23.9310 * (api / (degf + f2r)))
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
            pbmethod=pbmethod,
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
        # co = -1/bo*(dbodp - bg*drsdp/5.61458)
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
            gas_bg(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
            / 5.61458
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
            (0.00302 + 1.505 * rho_bs ** -0.951) * (degf - tscf) ** 0.938
            - (0.0216 - 0.0233 * 10 ** (-0.0161 * rho_bs))
            * (degf - tscf) ** 0.475
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
    pbmethod: pb_method = pb_method.VALMC,
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
    wt: Salt wt% (0-100). Default = 0
    ch4_sat: Degree of methane saturation (0 - 1). Default = 0
    export: Boolean flag that controls whether to export full table to excel, and separate PVDG and PVDO include files. Default is False
    pvto: Boolean flag that controls whether the pvto live oil Eclipse format will be generated. This can only be active if export flag is also True;
          - extends bubble point line up to maximum pressure
          - generates undersaturated oil propeties
          - writes out pvto include file
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
            sg_sp=sg_g,
            pb=pb,
            pbmethod=pbmethod,
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
            print(
                "Iteratively solving for Rsb fraction to use in order to harmonize user specified Pb and Rsb\n"
            )
            pbcalc = oil_pbub(
                degf=degf, api=api, sg_sp=sg_g, rsb=rsb, pbmethod=pbmethod
            )
            err = 100
            rsb_old = rsb
            i = 0
            while err > 0.0001:
                rsbnew = pb / pbcalc * rsb_old
                pbcalc = oil_pbub(
                    degf=degf,
                    api=api,
                    sg_sp=sg_g,
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
        print(
            "Iteratively solving for Rsb fraction to use at maximum pressure to deliver appropriate Pb and Rsb\n"
        )
        rsb_max = oil_rs_bub(
            degf=degf,
            api=api,
            sg_sp=sg_g,
            pb=pmax,
            pbmethod=pbmethod,
            rsmethod=rsmethod,
        )
        rs_at_pbi = oil_rs(
            api=api,
            degf=degf,
            sg_sp=sg_g,
            p=pb,
            pb=pmax,
            rsb=rsb_max,
            rsmethod=rsmethod,
            pbmethod=pbmethod,
        )
        err = rs_at_pbi - rsb
        rsb_old = rs_at_pbi
        rsb_frac_new = rsb_frac
        i = 0
        while abs(err) > 0.0001:
            rsb_frac_new = rsb / rs_at_pbi * rsb_frac
            rs_at_pbi = oil_rs(
                api=api,
                degf=degf,
                sg_sp=sg_g,
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
                    sg_sp=sg_g,
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
                sg_sp=sg_g,
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
                sg_g=sg_sp,
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
            gas_bg(p=p, sg=sg_sp, degf=degf, zmethod=zmethod, cmethod=cmethod)
            * 1000
            / 5.61458
        )  # rb/mscf
        gz.append(
            gas_z(p=p, sg=sg_sp, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        visg.append(
            gas_ug(p=p, sg=sg_sp, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        cg.append(gas_cg(p=p, sg=sg_sp, degf=degf, cmethod=cmethod))
        bw, lden, visw, cw, rsw = brine_props(
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
                            sg_sp=sg_g,
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

    st_deno = sg_o * 62.4  # lb/cuft
    st_deng = gas_den(
        p=psc, sg=sg_sp, degf=tscf, zmethod=zmethod, cmethod=cmethod
    )
    bw, lden, visw, cw, rsw = brine_props(
        p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat
    )
    res_denw = lden * 62.4  # lb/cuft
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
        df.to_excel("bot.xlsx", engine="openpyxl")
        pvdg = df[["Pressure (psia)", "Bg (rb/mscf", "ug (cP)"]]
        pvdg = pvdg.set_index("Pressure (psia)")
        headers = ["-- P (psia)", "Bg (rb/mscf", "ug (cP)"]
        fileout = "PVDG\n" + tabulate(pvdg, headers) + "\n/"
        with open("PVDG.INC", "w") as text_file:
            text_file.write(fileout)
        pvdo = df[["Pressure (psia)", "Bo (rb/stb)", "uo (cP)"]]
        pvdo = pvdo.set_index("Pressure (psia)")
        headers = ["-- P (psia)", "Bo (rb/stb)", "uo (cP)"]
        fileout = "PVDO\n" + tabulate(pvdg, headers) + "\n/"
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


def gas_water_content(p: float, degf: float) -> float:
    """ Returns saturated volume of water vapor in natural gas (stb/mmscf)
        From 'PVT and Phase Behaviour Of Petroleum Reservoir Fluids' by Ali Danesh
        degf: Water Temperature (deg F)
        p: Water pressure (psia)
    """
    t = degf
    content = (
        (
            47484
            * (
                np.exp(
                    69.103501
                    + (-13064.76 / (t + f2r))
                    + (-7.3037 * np.log(t + f2r))
                    + (0.0000012856 * ((t + f2r) * (t + f2r)))
                )
            )
            / (p)
            + (np.power(10, ((-3083.87 / (t + f2r)) + 6.69449)))
        )
        * (1 - (0.00492 * 0) - (0.00017672 * (0 * 0)))
        / 8.32
        / 42
    )
    return content


def brine_props(p: float, degf: float, wt: float, ch4_sat: float) -> tuple:
    """ Calculates Brine properties from modified Spivey Correlation per McCain Petroleum Reservoir Fluid Properties pg 160
        Returns tuple of (Bw (rb/stb), Density (sg), viscosity (cP), Compressibility (1/psi), Rw GOR (scf/stb))
        p: Pressure (psia)
        degf: Temperature (deg F)
        wt: Salt wt% (0-100)
        ch4_sat: Degree of methane saturation (0 - 1)
    """

    def Eq41(
        t, input_array
    ):  # From McCain Petroleum Reservoir Fluid Properties
        t2 = t / 100
        return (
            input_array[1] * t2 ** 2 + input_array[2] * t2 + input_array[3]
        ) / (input_array[4] * t2 ** 2 + input_array[5] * t2 + 1)

    Mpa = p * 0.00689476  # Pressure in mPa
    degc = (degf - 32) / 1.8  # Temperature in deg C
    degk = degc + 273  # Temperature in deg K
    m = (
        1000 * (wt / 100) / (58.4428 * (1 - (wt / 100)))
    )  # Molar concentration of NaCl from wt % in gram mol/kg water

    rhow_t70_arr = [0, -0.127213, 0.645486, 1.03265, -0.070291, 0.639589]
    Ewt_arr = [0, 4.221, -3.478, 6.221, 0.5182, -0.4405]
    Fwt_arr = [0, -11.403, 29.932, 27.952, 0.20684, 0.3768]
    Dm2t_arr = [0, -0.00011149, 0.000175105, -0.00043766, 0, 0]
    Dm32t_arr = [0, -0.0008878, -0.0001388, -0.00296318, 0, 0.51103]
    Dm1t_arr = [0, 0.0021466, 0.012427, 0.042648, -0.081009, 0.525417]
    Dm12t_arr = [0, 0.0002356, -0.0003636, -0.0002278, 0, 0]
    Emt_arr = [0, 0, 0, 0.1249, 0, 0]
    Fm32t_arr = [0, -0.617, -0.747, -0.4339, 0, 10.26]
    Fm1t_arr = [0, 0, 9.917, 5.1128, 0, 3.892]
    Fm12t_arr = [0, 0.0365, -0.0369, 0, 0, 0]

    rhow_t70 = Eq41(degc, rhow_t70_arr)
    Ewt = Eq41(degc, Ewt_arr)
    Fwt = Eq41(degc, Fwt_arr)
    Dm2t = Eq41(degc, Dm2t_arr)
    Dm32t = Eq41(degc, Dm32t_arr)
    Dm1t = Eq41(degc, Dm1t_arr)
    Dm12t = Eq41(degc, Dm12t_arr)
    Emt = Eq41(degc, Emt_arr)
    Fm32t = Eq41(degc, Fm32t_arr)
    Fm1t = Eq41(degc, Fm1t_arr)
    Fm12t = Eq41(degc, Fm12t_arr)

    cwtp = (1 / 70) * (1 / (Ewt * (Mpa / 70) + Fwt))  # Eq 4.2

    Iwt70 = (1 / Ewt) * np.log(abs(Ewt + Fwt))  # Eq 4.3
    Iwtp = (1 / Ewt) * np.log(abs(Ewt * (Mpa / 70) + Fwt))  # Eq 4.4
    rhowtp = rhow_t70 * np.exp(Iwtp - Iwt70)  # Eq 4.5

    rhobt70 = (
        rhow_t70
        + Dm2t * m * m
        + Dm32t * m ** 1.5
        + Dm1t * m
        + Dm12t * m ** 0.5
    )  # Eq 4.6
    Ebtm = Ewt + Emt * m  # Eq 4.7
    Fbtm = Fwt + Fm32t * m ** 1.5 + Fm1t * m + Fm12t * m ** 0.5  # Eq 4.8
    cbtpm = (1 / 70) * (1 / (Ebtm * (Mpa / 70) + Fbtm))  # Eq 4.9
    Ibt70 = (1 / Ebtm) * np.log(abs(Ebtm + Fbtm))  # Eq 4.10
    Ibtpm = (1 / Ebtm) * np.log(abs(Ebtm * (Mpa / 70) + Fbtm))  # Eq 4.11
    Rhob_tpm = rhobt70 * np.exp(
        Ibtpm - Ibt70
    )  # Eq 4.12 - Density of pure brine (no methane) in SG

    # Re-evaluate at standard conditions (15 deg C)
    rhow_sc70 = Eq41(15, rhow_t70_arr)
    Ew_sc = Eq41(15, Ewt_arr)
    Fw_sc = Eq41(15, Fwt_arr)
    Dm2_sc = Eq41(15, Dm2t_arr)
    Dm32_sc = Eq41(15, Dm32t_arr)
    Dm1_sc = Eq41(15, Dm1t_arr)
    Dm12_sc = Eq41(15, Dm12t_arr)
    Em_sc = Eq41(15, Emt_arr)
    Fm32_sc = Eq41(15, Fm32t_arr)
    Fm1_sc = Eq41(15, Fm1t_arr)
    Fm12_sc = Eq41(15, Fm12t_arr)

    cw_sc = (1 / 70) * (1 / (Ew_sc * (0.1013 / 70) + Fw_sc))
    Iw_sc70 = (1 / Ew_sc) * np.log(abs(Ew_sc + Fw_sc))
    Iw_sc = (1 / Ew_sc) * np.log(abs(Ew_sc * (0.1013 / 70) + Fw_sc))
    rhow_sc = rhow_sc70 * np.exp(Iw_sc - Iw_sc70)
    rhob_sc70 = (
        rhow_sc70
        + Dm2_sc * m * m
        + Dm32_sc * m ** 1.5
        + Dm1_sc * m
        + Dm12_sc * m ** 0.5
    )
    Eb_scm = Ew_sc + Em_sc * m
    Fb_scm = Fw_sc + Fm32_sc * m ** 1.5 + Fm1_sc * m + Fm12_sc * m ** 0.5
    cb_scm = (1 / 70) * (1 / (Eb_scm * (0.1015 / 70) + Fb_scm))
    Ib_sc70 = (1 / Eb_scm) * np.log(abs(Eb_scm + Fb_scm))
    Ib_scm = (1 / Eb_scm) * np.log(abs(Eb_scm * (0.1015 / 70) + Fb_scm))
    Rhob_scm = rhob_sc70 * np.exp(
        Ib_scm - Ib_sc70
    )  # Density of pure brine (no methane) in SG at standard conditions

    a_coefic = [
        0,
        -7.85951783,
        1.84408259,
        -11.7866497,
        22.6807411,
        -15.9618719,
        1.80122502,
    ]
    x = 1 - (degk / 647.096)  # Eq 4.14
    ln_vap_ratio = (647.096 / degk) * (
        a_coefic[1] * x
        + a_coefic[2] * x ** 1.5
        + a_coefic[3] * np.power(x, 3)
        + a_coefic[4] * np.power(x, 3.5)
        + a_coefic[5] * np.power(x, 4)
        + a_coefic[6] * np.power(x, 7.5)
    )  # Eq 4.13
    vap_pressure = np.exp(ln_vap_ratio) * 22.064

    a_coefic = [0, 0, -0.004462, -0.06763, 0, 0]
    b_coefic = [0, -0.03602, 0.18917, 0.97242, 0, 0]
    c_coefic = [0, 0.6855, -3.1992, -3.7968, 0.07711, 0.2229]

    A_t = Eq41(degc, a_coefic)
    B_t = Eq41(degc, b_coefic)
    C_t = Eq41(degc, c_coefic)

    try:
        mch4w = np.exp(A_t * np.power(np.log(Mpa - vap_pressure), 2) + B_t * np.log(Mpa - vap_pressure) + C_t)  # Eq 4.15
    except:
        mch4w = 0
    
    u_arr = [
        0,
        8.3143711,
        -7.2772168e-4,
        2.1489858e3,
        -1.4019672e-5,
        -6.6743449e5,
        7.698589e-2,
        -5.0253331e-5,
        -30.092013,
        4.8468502e3,
        0,
    ]
    lambda_arr = [
        0,
        -0.80898,
        1.0827e-3,
        183.85,
        0,
        0,
        3.924e-4,
        0,
        0,
        0,
        -1.97e-6,
    ]
    eta_arr = [0, -3.89e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    lambda_ch4Na = (
        lambda_arr[1]
        + lambda_arr[2] * degk
        + (lambda_arr[3] / degk)
        + lambda_arr[6] * Mpa
        + lambda_arr[10] * Mpa * Mpa
    )
    Eta_ch4Na = eta_arr[1]
    mch4b = mch4w * np.exp(
        -2 * lambda_ch4Na * m - Eta_ch4Na * m * m
    )  # Eq 4.18 - Methane solubility in brine (g-mol/kg H2O)

    mch4 = ch4_sat * mch4b  # Fraction of saturated methane solubility

    dudptm = (
        u_arr[6]
        + u_arr[7] * degk
        + (u_arr[8] / degk)
        + (u_arr[9] / (degk * degk))
    )  # Eq 4.19
    dlambdadptm = lambda_arr[6] + 2 * lambda_arr[10] * Mpa  # Eq 4.20
    detadptm = 0  # Eq 4.21

    Vmch4b = (
        8.314467 * degk * (dudptm + 2 * m * dlambdadptm + m * m * 0)
    )  # Eq 4.22
    vb0 = 1 / Rhob_tpm  # Eq 4.23
    rhobtpbch4 = (1000 + m * 58.4428 + mch4 * 16.043) / (
        (1000 + m * 58.4428) * vb0 + (mch4 * Vmch4b)
    )  # Eq 4.24... mch4 = Methane concentration in g/cm3
    vbtpbch4 = 1 / rhobtpbch4
    dvbdp = -vb0 * cbtpm  # Eq 4.27
    d2uch2dp2 = 0
    d2lambdadp2 = 2 * lambda_arr[10]
    d2etadp2 = 0
    dVmch4dp = (
        8.314467 * degk * (d2uch2dp2 + 2 * m * d2lambdadp2 + m * m * d2etadp2)
    )  # Eq 4.31
    cwu = -((1000 + m * 58.4428) * dvbdp + mch4 * dVmch4dp) / (
        (1000 + m * 58.4428) * vb0 + (mch4 * Vmch4b)
    )  # Eq 4.32 -- Undersaturated brine Compressibility (Mpa-1)
    satdmch4dp = (
        mch4
        * (2 * A_t * np.log(Mpa - vap_pressure) + B_t)
        / ((Mpa - vap_pressure) - 2 * dlambdadptm * m)
    )  # Eq 4.33

    zee = gas_z(p=p, sg=0.5537, degf=degf)  # Z-Factor of pure methane

    vmch4g = zee * 8.314467 * degk / Mpa  # Eq 4.34

    cws = -(
        (1000 + m * 58.4428) * dvbdp
        + mch4 * dVmch4dp
        + satdmch4dp * (Vmch4b - vmch4g)
    ) / (
        (1000 + m * 58.4428) * vb0 + (mch4 * Vmch4b)
    )  # Eq 4.35 - Compressibility of saturated brine Mpa-1
    cw_new = 1 / (145.038 * (1 / cws))  # Compressibility in psi-1
    vb0_sc = (
        1 / Rhob_scm
    )  # vb0 at standard conditions - (Calculated by evaluating vbo at 0.1013 MPa and 15 degC)
    Bw = (((1000 + m * 58.4428) * vb0) + (mch4 * Vmch4b)) / (
        (1000 + m * 58.4428) * vb0_sc
    )

    # m =  Molar concentration of NaCl from wt % in gram mol/kg water
    # mch4b = Methane solubility in brine (g-mol/kg H2O)
    # mch4 = Fraction of saturated methane solubility
    # Vmch4b = 



    zee_sc = gas_z(p=psc, sg=0.5537, degf=tscf)
    vmch4g_sc = zee_sc * 8.314467 * (273 + 15) / 0.1013  # Eq 4.34
    rsw_new = mch4 * vmch4g_sc / ((1000 + m * 58.4428) * vb0_sc)
    rsw_new_oilfield = rsw_new / 0.1781076  # Convert to scf/stb

    d = [
        0,
        2885310,
        -11072.577,
        -9.0834095,
        0.030925651,
        -0.0000274071,
        -1928385.1,
        5621.6046,
        13.82725,
        -0.047609523,
        0.000035545041,
    ]
    a = [-0.21319213, 0.0013651589, -0.0000012191756]
    b = [0.069161945, -0.00027292263, 0.0000002085244]
    c = [-0.0025988855, 0.0000077989227]

    lnuw_tp = sum([d[i] * np.power(degk, (i - 3)) for i in range(1, 6)])
    lnuw_tp += sum(
        [rhowtp * (d[i] * np.power(degk, (i - 8))) for i in range(6, 11)]
    )

    uw_tp = np.exp(lnuw_tp)

    AA = a[0] + a[1] * degk + a[2] * degk * degk  # Eq 4.43
    BB = b[0] + b[1] * degk + b[2] * degk * degk
    CC = c[0] + c[1] * degk

    lnur_tm = AA * m + BB * m * m + CC * m * m * m  # Eq 4.46
    ur_tm = np.exp(lnur_tm)
    ub_tpm = ur_tm * uw_tp * 1000  # cP - Eq 4.48

    bw = Bw  # rb/stb
    lden = rhobtpbch4  # sg
    visw = ub_tpm  # cP
    cw = cw_new  # 1/psi
    rsw = rsw_new_oilfield  # scf/stb

    return (bw, lden, visw, cw, rsw)


def lorenz2b(lorenz: float, lrnz_method: str = "EXP") -> float:
    """ Returns B-factor that characterizes the Lorenz function
        Lorenz: Lorenz coefficient (0-1)
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' (Exponential) or 'LANG' (Langmuir).
                Defaults to EXP if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
    """
    method = lrnz_method.upper()
    if method != "EXP" and method != "LANG":
        print('Method must be "LANG" or "EXP"')
        sys.exit()

    if lorenz < 0.000333:
        B = 2 / 1000
        if method == "LANG":
            B = 1 / 1000
        return B
    if lorenz > 0.997179125528914:
        B = 709
        if method == "LANG":
            B = 25000
        return B

    # Set bookends for B
    hi = 709
    if method == "LANG":
        hi = 25000
    lo = 0.000001
    args = (lorenz, method)

    def LorenzErr(args, B):
        lorenz, method = args
        B = max(B, 0.000001)
        if method == "EXP":
            B = min(B, 709)
            err = 2 * ((1 / (np.exp(B) - 1)) - (1 / B)) + 1 - lorenz
        else:
            B = min(B, 25000)
            PL = 1 / B
            VL = PL + 1
            err = (
                VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5
            ) * 2 - lorenz
        return err

    rtol = 0.0000001
    return bisect_solve(args, LorenzErr, lo, hi, rtol)


def lorenzfromb(B: float, lrnz_method: str = "EXP") -> float:
    """ Returns Lorenz coefficient that corresponds to a Beta value
        B: The B-Factor (positive float)
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
    """
    method = lrnz_method.upper()
    B = max(B, 0.000001)
    if method == "LANG":
        B = min(B, 25000)
        PL = 1 / B
        VL = PL + 1
        L = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
    else:
        B = min(B, 709)
        L = 2 * (1 / (np.exp(B) - 1) - (1 / B)) + 1
    return L


def lorenz_from_flow_fraction(
    kh_frac: float, phih_frac: float, lrnz_method: str = "EXP"
) -> float:
    """ Returns Lorenz coefficient consistent with observed flow fraction from a phi_h fraction
        kh_frac: (0 - 1). Fraction of total flow from best quality reservoir phi_h
        phih_frac: (0 - 1). phi_h fraction that delivers the observed kh_fraction of flow
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
    """
    method = lrnz_method.upper()
    if kh_frac <= phih_frac:  #
        print("kh fraction should always be greater than phi_h fraction")
        return 0.001
    if kh_frac >= 1:
        print("kh Fraction must be less than 1")
        return 0.001

    # If Langmuir method, can explicitly calculate B
    if method == "LANG":
        x = phih_frac
        y = kh_frac
        B = (y - x) / (x * (1 - y))
        return lorenzfromb(B, method)

    # Set bookends and first guess of B
    hi = 709
    lo = 0.000001
    args = (kh_frac, phih_frac, method)

    def BErr(args, B):
        kh_frac, phih_frac, method = args
        method = method.upper()
        B = max(B, 0.000001)
        if method == "EXP":
            B = min(B, 709)
            err = (1 - np.exp(-B * phih_frac)) / (1 - np.exp(-B)) - kh_frac
        else:
            B = min(B, 25000)
            PL = 1 / B
            VL = PL + 1
            err = (VL * phih_frac) / (PL + phih_frac) - kh_frac
        return err

    rtol = 0.0000001
    B = bisect_solve(args, BErr, lo, hi, rtol)
    return lorenzfromb(B, method)


def lorenz_2_flow_frac(
    lorenz: float, phih_frac: float, lrnz_method: str = "EXP", B: float = -1
) -> float:
    """ Returns expected flow fraction from the best phi_h fraction, with a specified Lorenz coefficient

        lorenz: (0-1) Lorenz hetrogeneity factor
        phih_frac: (0 - 1). Best phi_h fraction
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
        B: Factor that characterizes the Lorenz function for the given method. Will calculate if only lorenz variable defined
        lorenz: Lorenz coefficient (0-1). If B is provided, will ignore this parameter to be more efficient. If not, will calculate B from this parameter.
    """

    method = lrnz_method.upper()
    if B < 0 and lorenz < 0:
        print("Must define either B or lorenz parameters")
        sys.exit()

    if B < 0:  # Need to calculate B
        B = lorenz2b(lorenz=lorenz, lrnz_method=lrnz_method)

    B = max(B, 0.000001)
    if method == "EXP":
        B = min(B, 709)
        fraction = (1 - np.exp(-B * phih_frac)) / (1 - np.exp(-B))
    else:
        B = min(B, 25000)
        PL = 1 / B
        VL = PL + 1
        fraction = (VL * phih_frac) / (PL + phih_frac)
    return fraction


def lorenz_2_layers(
    lorenz: float,
    k_avg: float,
    nlayers: int = 1,
    shuffle: bool = False,
    lrnz_method: str = "EXP",
    B: float = -1,
    phi_h_fracs: list = [],
) -> np.ndarray:
    """ Returns np.array of permeability values honoring a specified average permeability (assuming equal thickness layers unless list of phi_h_fracs is provided), with degree of heterogeneity consistant with specified Lorenz coefficient and method

        If B is left default, then it will be calculated. If B is explictly specified > 0, then it will be used instead of the provided lorenz coefficient so as to eliminate repetitive solving for B.


        lorenz: Lorenz coefficient (0-1). If B is provided, will igonore this parameter to be more efficient. If not, will calculate B from this parameter.
        nlayers: The number of permeability layers desired (>1 needed unless a list of phi_h_fracs is supplied)
        kavg: The average permeability of all the layers (assuming equal thickness)
        shuffle: Boolean flag to determine whether to return the permeability array in decreasing order (False), or random order (True). Default False. Will be reset to False if user defined phi_h_fracs are supplied
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
        B: Factor that characterizes the Lorenz function for the given method. Will calculate if only lorenz variable defined
        phi_h_fracs: Optional ability to specify a sorted list of phi_h fractions to get permeabilities for. If this list does not add to unity, then one additional layer permeability will be returned. The list needs to be in sorted order of best flow capacity to worst

    """
    if nlayers <= 1:
        if len(phi_h_fracs) < 2:
            return np.array([k_avg])

    method = lrnz_method.upper()

    if B < 0:  # Need to calculate B
        B = lorenz2b(lorenz=lorenz, lrnz_method=lrnz_method)

    B = max(B, 0.000001)
    if method == "EXP":
        B = min(B, 709)
    else:
        B = min(B, 25000)

    user_layers = False
    if len(phi_h_fracs) > 1:
        user_layers = True
        if sum(phi_h_fracs) > 1:
            phi_h_fracs = [x / sum(phi_h_fracs) for x in phi_h_fracs]
        if sum(phi_h_fracs) < 1:
            phi_h_fracs.append(1 - sum(phi_h_fracs))
        phih = (
            [0]
            + [sum(phi_h_fracs[: i + 1]) for i in range(len(phi_h_fracs) - 1)]
            + [1.0]
        )
        nlayers = len(phi_h_fracs)
    else:
        phih = np.arange(0, 1 + 1 / (nlayers), 1 / (nlayers))
        phi_h_fracs = np.array([1 / nlayers for i in range(len(phih) - 1)])
    sumkh = []

    for layer in phih:
        if method == "EXP":
            sumkh.append((1 - np.exp(-B * layer)) / (1 - np.exp(-B)))
        else:
            PL = 1 / B
            VL = PL + 1
            sumkh.append((VL * layer) / (PL + layer))

    kh = (
        np.array([sumkh[i] - sumkh[i - 1] for i in range(1, len(sumkh))])
        * k_avg
    )
    k = kh / np.array(phi_h_fracs)
    if shuffle:
        if not user_layers:
            np.random.shuffle(k)
    return k


class component_library:
    def __init__(self, model='PR79'):
        path = 'component_library.xlsx'
        filepath = pkg_resources.resource_filename(__name__, path)
        self.df = pd.read_excel(filepath, engine="openpyxl")
        self.model = model
        self.all_cols = ['Name', 'MW', 'Tc_R', 'Pc_psia',
                         'Visc_Zc', 'Pchor', 'Vc_cuft_per_lbmol']
        self.model_cols = ['Acentric', 'VTran', 'Tb_F', 'SpGr']
        self.all_dics = {}
        self.model_dics = {}
        self.components = self.df['Component'].tolist()
        self.names = self.df['Name'].tolist()
        self.property_list = self.all_cols + self.model_cols
        # Create dictionaries for all the model agnostic properties
        for col in self.all_cols:
            self.all_dics[col.upper()] = dict(
                zip(self.df['Component'], self.df[col]))
        # And then for all the model specific properties
        self.models = ['PR79', 'PR77', 'SRK', 'RK']
        for model in self.models:
            model_dic = {}
            for col in self.model_cols:
                model_dic[col.upper()] = dict(
                    zip(self.df['Component'], self.df[model+'-'+col]))
            self.model_dics[model] = model_dic

    def prop(self, comp, prop, model='PR79'):
        comp = comp.upper()
        if comp not in self.components:
            return 'Component not in Library'
        prop = prop.upper()
        props = [x.upper() for x in self.property_list]
        if model.upper() not in self.models:
            return 'Incorrect Model Name'
        dic = self.model_dics[model.upper()]
        if prop == 'ALL':
            return [self.all_dics[p.upper()][comp] for p in self.all_cols]+[dic[p.upper()][comp] for p in self.model_cols]
        props = [x.upper() for x in self.property_list]
        if prop not in props:
            return 'Property not in Library'
        if prop in [x.upper() for x in self.all_cols]:
            return self.all_dics[prop][comp]
        if prop in [x.upper() for x in self.model_cols]:
            return dic[prop][comp]
        return 'Component or Property not in library'

comp_library = component_library()

# CO2:Brine Library coded by Mark Burgoyne, July 2023
# Leveraging some translated SRK code snippets from VBA supplied by Steve Furnival

#import numpy as np
#import math
#from matplotlib import pyplot as plt

#import warnings
#warnings.filterwarnings("error")

EPS = 1e-8

#--Mole Weights ----------------------
MWSAL = 58.4428        # Mole Weight of Salt (NaCl)
MWWAT = 18.01528       # Mole Weight of Pure Water
MWCO2 = 44.01          # Mole Weight of CO2

#--Unit Conversions & Constants--------------------------------------------
BAR2PSI = 14.5037738
RGASCON = 83.1447      # Units of bar.cm3/(mol.K)
CEL2KEL = 273.15       # deg K at 0 degC
CONMOLA = 1000 / MWWAT # Moles in 1000kg Water [Molality Conversion Factor]
PSTND = 1.01325        # Standard Pressure (bar)
TSTND = 288.706        # Standard Temperature (Kelvin)
VMCO2S = 23690.5       # Molar Volume of CO2 at STP [cm3/gmol]
RHOCO2S = 0.00185771   # Density CO2 at STP [gm/cm3]
DENW = 998.98          # Freshwater density at standard conditions [kg/m3]
KGMOL2SM3 = 23.545     # sm3 CO2 per kg-mol at 60 deg F and 1 atm (from PhazeComp run to get definitive value, where zCO2 = 0.99388)
BBL2CUFT = 5.614583333 # cuft in a bbl

#============================================================================
#    ***  Mutual solubilities between CO2 and Brine - Calculated with  ***
#  A Phase-Partitioning Model for CO2–Brine Mixtures at Elevated Temperatures 
#  and Pressures: Application to CO2-Enhanced Geothermal Systems
#  Nicolas Spycher & Karsten Pruess, Transp Porous Med (2010) 82:173–196
#  DOI 10.1007/s11242-009-9425-y
#============================================================================

#===================================================================================================================
#  ***  Base Brine properties calculated using Spycher calculated xCO2. Additional calculations via the following ***
#Property                       Source
#-----------------------------  --------------------------------------------------------------------------------------------
#Pure Brine Density             Spivey et al. (modified),
#                               per "Petroleum Reservoir Fluid Property Correlations", (McCain, Spivey & Lenn: Chapter 4)
#CO2 Corrected Brine Density    Molar volume of dissolved CO2 estimated with Garcia (2001) equation, used with xCO2 calculated 
#                               from Spycher & Pruess, and CO2-free brine density from Spivey et al to calculate insitu density
#Pure Brine viscosity           Mao-Duan (2009) approach for pure brine viscosity
#CO2 Corrected Brine Viscosity  Used approach from "Viscosity Models and Effects of Dissolved CO2", Islam-Carlson (2012)
#                               to adjust the pure brine viscosity for xCO2 calculated from Spycher & Pruess
#===================================================================================================================

class CO2_Brine_Mixture():
    """ Calculates CO2 saturated Brine mutual solubilities and brine properties
    
            Inputs:
                pres: Pressure (Bar / psia)
                temp: Temperature (deg C / deg F)
                ppm: NaCL equivalent weight concentration in brine in parts NaCl per million parts of brine (default zero, Wt% = 100 * ppm / 1E6 )
                metric: Boolean operator that determines units assumed for input, and return calculated (default True)
                cw_sat: Boolean operator that determines whether to calculate saturated brine compressibility, doubling calculations required (default False)
    
            Returns object with following calculated properties:
                .x       : Mole fractions of CO2 and H2O in aqueous phase [xCO2, xH2O]
                .y       : Mole fractions of CO2 and H2O in vapor phase [yCO2, yH2O]
                .xSalt   : Sum of Mole fraction of Na + Cl species in brine (Double the NaCl species mole fraction)
                .rhoGas  : CO2 rich gas density (gm/cm3)
                .bDen    : Brine density (gm/cm3) [CO2 Saturated, Pure Brine, Freshwater]
                .bVis    : Brine viscosity (cP)   [CO2 Saturated, Pure Brine, Freshwater]
                .bVisblty: CO2 Saturated brine viscosibility (1/Bar or 1/psi)
                .bw      : Brine formation volume factor (rm3/smr / rb/stb) [CO2 Saturated, Pure Brine, Freshwater]
                .Rs      : CO2 Saturated Brine solution gas ratio (sm3/sm3 or scf/stb), relative to standard conditions
                .Cf_usat : Brine undersaturated compressibility (1/Bar or 1/psi). The compressibility with constant Rs
                .Cf_sat  : Brine saturated compressibility (1/Bar or 1/psi). The compressibility with reducing Rs under depletion
                            
            Usage example for 5000 psia x 275 deg F and 3% NaCl brine:
                mix = rtb.CO2_Brine_Mixture(pres = 5000, temp = 275, ppm = 30000, metric = False)
                mix.bw  # Returns [CO2 Saturated Brine Bw, Pure Brine Bw, Pure Water Bw]
                >> [1.1085795290443725, 1.0543051245909865, 1.0542061001251017]
                
                mix.x  # Returns [xCO2, xBrine]
                >> array([0.02431245, 0.95743175])
                
            Usage example for 175 Bara x 85 degC and 0% NaCl brine:
                mix = rtb.CO2_Brine_Mixture(pres = 175, temp = 85)
                mix.Rs  # Returns sm3 dissolved CO2 / sm3 Brine
                >> 24.742923469934272
                
           
    """
    def __init__(self, pres, temp, ppm = 0, metric = True, cw_sat = False):
        self.metric = metric              # Units. FIELD or METRIC
        self.ppm = ppm                    # Parts (by wt) NaCl added to 1E6 parts of water
        self.tKel = None                  # Deg K
        self.P0 = None                    # Standard pressure (Bar)
        self.fugPi = [None, None]         # Fugacity * Pressure for species i
        self.K = [None, None]             # yi/xi at reservoir pressure
        self.gamma_prime = None
        self.gamma = ([1, 1])
        self.pRT = None                   # P/RT
        self.pRT0 = None                  # (P - P0)/RT
        self.y = np.array([1.0, 0.0])     # Mole fraction split CO2 and H2O in vapor phase
        self.x = np.array([0.0, 1.0])     # Mole fraction split CO2 and H2O in aqueous phase
        self.A = None
        self.Bprime = None
        self.xSalt = None                 # Mole fraction salt in aqueous mixture
        self.molaL = None                 # Brine Molality (gmol/kg)
        self.MolarVol = None              # Molar Volumes from RK-EOS (cm3/gmol)
        self.GASZ = None                  # Vapor phase Z-Factor
        self.MwGas = None                 # Mole weight  [gm/gmol]
        self.MwBrine = None               # Mole weight of CO2 free brine (gm/gmol)
        self.rhoGas = None                # Density of Gas Mixture (gm/cm3)
        self.aMix = None                  # RK-EOS A_mix parameter
        self.aij = None                   # RK-EOS aij
        self.kij = None                   # RK-EOS kij
        self.bMix = None                  # RK-EOS B_mix parameter
        self.b = None                     # RK-EOS bij
        self.vBar = [0, 0]                # Avg partial molar volume of pure condensed phase over the pressure interval P0 to pBar
        self.bDen = None                  # CO2 laden Brine density (gm/cm3) [CO2 Saturated, CO2 Free]
        self.bVis = None                  # Brine viscosity (cP)  [CO2 Saturated, CO2 Free]
        self.bVisblty = None              # CO2 laden viscosibility (1/Bar or 1/psi) (at pressures above Psat)
        self.bw = None                    # Brine formation volume factor (res vol/std vol) [CO2 Saturated, CO2 Free]
        self.Rs = None                    # Solution gwr CO2 (sm3/sm3 or scf/stb brine) relative to standard conditions
        self.Cf_usat = None               # Undersaturated brine compressibility (1/Bar or 1/psi). Compressibility without changing Rs
        self.Cf_sat = None                # Saturated brine compressibility (1/Bar or 1/psi). Compressibility with changing Rs
        self.CO2_sat = False              # Flag to determine if in P-T range for saturated liquid CO2 K values
        self.repeat = False               # Flag to trigger repeat of calculations depending on whether CO2 is saturated liquid phase
        #self.EzrokhiDenA = None           # Ezrokhi coefficient array for density (to be calculated with ezrokhi() function)
        #self.EzrokhiVisB = None           # Ezrokhi coefficient array for viscosity (to be calculated with ezrokhi() function)
        self.Rs_STD = None                # Dissolved CO2 remaining at standard conditions (sm3/sm3)
        self.ppm_sat = None               # Maximum ppm Salt at specified temperature
        
        xNaCl = (ppm / MWSAL) / ((ppm / MWSAL) + (1000000 - ppm) / MWWAT) # Mole fraction of salt in pure brine (note, this is not the same as self.xSalt used in EOS calculations as those are actually xNa + xCl, ie double the NaCl species mole fraction)
        self.MwBrine = xNaCl * MWSAL + (1 - xNaCl) * MWWAT
        
        if self.metric:
            self.pBar = pres              # Pressure (Bar)
            self.degC = temp              # Deg C
        else:
            self.pBar = pres / BAR2PSI      # Pressure (psia -> Bar)
            self.degC = (temp - 32)/1.8   # Temperature (degF -> deg C)
        self.tKel = self.degC + CEL2KEL
        
        # Calculate maximum salt concentration
        self.ppm_sat = round(262180 + 72 * self.degC + 1.06 * self.degC**2,0)  # Eq 9.1 from Whitson Phase Monograph
        
        #if self.degC <= 31:
        #    # Below saturated CO2 pressure equation emprically fit to data at https://www.ohio.edu/mechanical/thermo/property_tables/CO2/CO2_TempSat1.html
        #    if self.pBar >= 10**-15.90106469 * self.tKel**7.157992919: # Above CO2 Psat 
        #        self.CO2_sat = True
        
        self.std_conditions_rs()          # Estimate Rs remain at standard conditions
        
        if cw_sat: # Determine properties at 0.5 Bar less pressure first, so that we can calculate saturated compressibility
            dP = 0.5
            self.pBar += dP
            self.co2BrineSolubility()         # Calculate mutual solubilities
            self.brine()                      # Calculate Brine properties 
            
            # Will ignore water in the vapor phase for compressibility calcs
            bw1 = self.bw[0]
            Rs1 = self.Rs # Ignore residual gas saturation, since we are looking at deltas, so they will cancel
            self.pBar -= dP

        self.co2BrineSolubility()         # Calculate mutual solubilities
        self.brine()                      # Calculate Brine properties       
        
        if cw_sat:
            bw2 = self.bw[0]
            Rs2 = self.Rs
            Bg = PSTND * self.GASZ * self.tKel / (TSTND * self.pBar)  # rm3/sm3 or rcf/scf
            
            dBwdP = (bw1 - bw2) / dP
            dRsdP = (Rs1 - Rs2) / dP
            if not self.metric:
                dRsdP /= BBL2CUFT  # Convert scf/stb to scf/scf = rm3/rm3
                
            self.Cf_sat = (1 / bw2)  *(-dBwdP + dRsdP * Bg)
            
            if not self.metric:
                self.Cf_sat /= BAR2PSI

    
    # Estimate residual Rs at standard conditions so that we dont have to repeat all Spycher calculations 
    # Rs at 1.01325 Bar and 288.706K was calculated over the ppm range of 0 - 290,000, and fit to the following
    # Equation by M. Burgoyne (Oct 2023), to offset the Spycher Rs to deliver zero Rs at standard conditions (equivalent to lab measurements)
    def std_conditions_rs(self):
        A, B, C, D = 0.1084E+03, -0.4573E+01, 0.7247E-06, -0.1433E+00
        # Equation of form: Y = A*EXP(B*EXP(C*X))+D <--- Gompertz + c
        self.Rs_STD = A * np.exp(B * np.exp(C * self.ppm)) + D
    
    # Ezrokhi functionality removed due to (a) not needed with full brine definitions and (b) some concern  as to how
    # IX is actually implementing subsequent density etc calculations (const CO2 molar volume which is incorrect).
    # Accordingly I've removed the functionality below in this script
    
        # Additional function available to calculate Ezrokhi coefficients for effects of dissolved CO2 on brine density and viscosity
        # .ezrokhi(lower_degC, upper_degC)
        #
        # Calculates and populates the following additional attributes;
        # .EzrokhiDenA : List of A_CO2's for equation Ai(T) = A[0] + A[1] * degC + A[2] * degC**2, for Ezrokhi density adjustment
        # .EzrokhiVisB : List of B_CO2's for equation Bi(T) = B[0] + B[1] * degC + B[2] * degC**2, for Ezrokhi viscosity adjustment
        
        #def ezrokhi(self, lower_degC, upper_degC):
        #    # Function to regress on Ezrokhi coefficients for the mixture over the temperature range
        #    # Will return arrays Ai and Bi for density and viscosity impact of dissolved CO2 in brine respectively
        #    # Uses given sample pressure
        #    
        #    # Ensure reversed temperatures, or cases with same upper & lower values don't fall over
        #    if lower_degC <= upper_degC: 
        #        lower_degC, upper_degC = min(lower_degC, upper_degC), max(lower_degC, upper_degC)
        #        lower_degC -= 10
        #        lower_degC = max(10, lower_degC) # Ensure doesn't get below 10 deg C
        #        upper_degC += 10
        #    
        #    
        #    temps = np.linspace(lower_degC, upper_degC)
        #    As, Bs = [], []
        #    
        #    # Grab the original temperature before doing calculations over range for Ezrokhi
        #    original_degc = self.degC 
        #    
        #    for temp in temps:
        #        self.degC = temp
        #        self.co2BrineSolubility()
        #        wt_co2 = self.x[0] * MWCO2 / (self.x[0] * MWCO2 + self.x[1] * self.MwBrine) # Weight fraction of dissolved CO2
        #        results = brine_props(self.pBar, temp, self.ppm, self.x[0], self.MwBrine)
        #        bDen, bVis, bVisblty, bw, Rs, Cf_usat = results
        #        As.append(np.log10(bDen[0]/bDen[1])/wt_co2)  # Corrections relative to CO2 free brine
        #        Bs.append(np.log10(bVis[0]/bVis[1])/wt_co2)  # Corrections relative to CO2 free brine
        #    
        #    # Because of the form of the Akand W. Islam and Eric S. Carlson viscosity adjustment, 
        #    # the Ezrokhi B coefficient for viscosity will be a constant value for a given xCO2
        #    self.EzrokhiVisB = [Bs[0], 0, 0]
        #    
        #    # Fit density data to quadratic form
        #    z = np.polyfit(temps, np.array(As), 2)
        #    # Numpy fits coefficients in reverse order to what is expected for Ezrokhi equation as follows
        #    # Ai = ao + a1 * T + a2 * T**2
        #    z = list(z)[::-1]  # Reverse order of fitted coefficients
        #    self.EzrokhiDenA = z 
        #    
        #    # And reset to original temperature & recalculate
        #    self.degC = original_degc
        #    self.co2BrineSolubility()
    
    def calc_type(self):
        # == Figure out which calculation method to employ ==========================
        # deg C <= 99: Non-iterative solution employed assuming negligible xCO2 for mixing rules
        # deg C > 109: An iterative solution is employed with more robust mixing rules application
        # 99 < deg C < 109: An iterative solution is employed. Partitioning factors and fugacity coefficients are blended between low temp and high temp relationships.
        self.low_temp = True      # Use the simple approach for temp < 99 deg C
        self.scaled = False       # No blending
        if self.degC > 99 and self.degC <= 109:
            self.low_temp = False # Use the iteretaive robust approach
            self.scaled = True    # Also use the simpler < 99 degC approach for equilibrium factors and fugacities and scale the result
        elif self.degC > 109:
            self.low_temp = False # Only perform the more robust > 99 degC approach
            self.scaled = False

    # Vapor pressure of water with Buck equation
    # https://en.wikipedia.org/wiki/Vapour_pressure_of_water#:~:text=The%20saturation%20vapour%20pressure%20of,pressure%20equals%20the%20ambient%20pressure.
    def water_vap_p(self):
        kpa = 0.61121 * np.exp((18.678 - (self.degC/234.5))*(self.degC/(257.14+self.degC)))
        return 0.145038 * kpa # psia
    
    # Initial estimate for yH2O in saturated Brine / CO2 system
    # Empirical fit by Mark Burgoyne against data generated with modified Whitson PR-EOS method, July 2023
    def est_yH2O(self):
        # Slope and intercept of ppH20/PPatm - 1
        # Slope as function of deg F
        # N.  98:   Y = (A+B*X)**C+D <--- Bleasdale (Shifted power) + c
        A, B, C, D = 0.1939E+01, 0.8913E-02, -.4844E+01, 0.2551E-03
        degf = self.degC * 1.8 + 32
        p = self.pBar * BAR2PSI
        slope = (A+B*degf)**C+D
        
        # Intercept as a function of degF
        # N. 107:   Y = -1/(A+B*X**2)**C <--- Inv pow
        A, B, C = 0.3097E+00, 0.9136E-05, 0.1719E+01
        intcpt = -1/(A+B*degf**2)**C
        
        pp_ratio = slope*p + intcpt + 1
        atm_pvap = self.water_vap_p()
        pp = pp_ratio * atm_pvap
        if pp < atm_pvap:
            pp = atm_pvap
        return max(min(pp/p, 1-EPS),0)
        
    def ppm2Molality(self):
    #=======================================================================
    #  Molality = gMoles Salt per kg of water 
    #  ppm / MWSAL = kgMoles Salt per 1,000,000 kg brine
    #  Subtract the ppm of salt to yield mass of water in 1,000,000 kg of brine in denominator
    #=======================================================================
        return self.ppm / MWSAL * 1000 / (1e6 - self.ppm)
    
    def blended_val(self, low_val, high_val): # Blends results between 99 - 109 deg C
        return ((self.degC - 99) * low_val + (109 - self.degC) * high_val)/(109 - 99)
        
    def aCO2_RK(self):
    #=======================================================================
    #  a-Coefficient of CO2 for RK-EoS is Temperature Dependent
    #  Delineating with low_temp flag to permit calculating low_temp relationship 
    #  values between 99 - 109 deg C
    #=======================================================================
        if self.low_temp:
            return self.FT(self.tKel, [7.54e7, -4.13e4])   # Low temperature relationship
        return self.FT(self.tKel, [8.008e7, -4.984e4])     # High temperature relationship
    
    def aH2O_RK(self):
    #=======================================================================
    # a-Coefficient of H2O for RK-EoS is Temperature Dependent
    # Only used for high temperature method
    #=======================================================================
        return self.FT(self.tKel, [1.337e8, -1.4e4])
    
    def gammaCO2(self):
        # Equation 18 from Spycher 2010 paper (Gamma dash)
    #=======================================================================
    #  CO2 Activity Coefficients
    #  Spycher, N., and Pruess, K.,
    #  "A Phase-Partitioning Model for CO2-Brine Mixtures ..."
    #  Transp Porous Med (2010), 82, pp. 173-196
    #=======================================================================
    
        cL = [0.0002217, 1.074, 2648.0]
        cZ = [0.000013, -20.12, 5259.0]
        
        #--(lamB,zetA) from Spycher & Pruess Equation (19) and Table 1---------------
        lamB = cL[0] * self.tKel + cL[1] / self.tKel + cL[2] / self.tKel ** 2
        zetA = cZ[0] * self.tKel + cZ[1] / self.tKel + cZ[2] / self.tKel ** 2
        
        #==Hassanzadeh Equation (12)============================================
        self.gamma_prime = (1.0 + self.molaL / CONMOLA) * np.exp(self.molaL * (2.0 * lamB + self.molaL * zetA))
    
    
    def aMix_RK(self):
    #=======================================================================
    #  a-Coefficient of CO2-H2O Mixture for RK-EoS
    #=======================================================================
        
        # Setup Kij array    
        K = np.zeros((2,2))
        if not self.low_temp:  # Use high temp method
            K[0][1] = self.FT(self.tKel, [0.4228, -7.422e-4])   # KCO2-H2O
            K[1][0] = self.FT(self.tKel, [1.427e-2, -4.037e-4]) # KH2O-CO2
        else:
            K[0][1] = 7.89e7
            K[1][0] = 7.89e7
        
        # Calculate kij array
        k = np.zeros((2,2))
        k[0][1] = (K[0][1] * self.y[0]) + (K[1][0] * self.y[1])  # Eq. A-6
        k[1][0] = k[0][1]
        
        # Setup aij array
        a = np.zeros((2,2))
        a[0][0] = self.aCO2_RK()
        a[1][1] = self.aH2O_RK()
        if self.low_temp:
            a[0][1] = 7.89e7
            a[1][0] = a[0][1]
        else:
            for i in range(2):
                j = 1 - i
                a[i][j] = (a[i][i] * a[j][j])**0.5 *(1-k[i][j])  # Eq. A-5
        amix = 0
        for i in range(2):
            for j in range(2):
                amix += self.y[i]*self.y[j]*a[i][j]
            
        self.aMix = amix
        self.aij = a
        self.kij = k
        
    
    def bMix_RK(self):
    #=======================================================================
    #  b-Coefficient of CO2-H2O Mixture for RK-EoS
    #=======================================================================
        if self.low_temp:
            b = np.array([27.80, 18.18])
        else:
            b = np.array([28.25, 15.70])
        
        self.bMix = np.dot(self.y, b)
        self.b = b
    
    def cubicSolver(self, e2, e1, e0):
    #=======================================================================
    #  Cubic Polynomial Solver: f(Z) = Z**3 + E2*Z**2 + E1*Z + E1 = 0
    #=======================================================================
        self.repeat = False
        Z = np.roots(np.array([1.0, e2, e1, e0]))
        Z = np.array([x for x in Z if np.isreal(x)]) # Keep only real results
        if len(Z) > 1: # Evaluate which root to use per Eqs 25 and 26 in Spycher & Pruess (2003)
            vgas, vliq = max(Z), min(Z)
            
            w1 = self.pBar*(vgas - vliq)
            w2 = RGASCON * self.tKel * np.log((vgas - self.bMix)/(vliq - self.bMix)) + self.aMix/(self.tKel**0.5 * self.bMix) * np.log((vgas + self.bMix) * vliq / ((vliq + self.bMix) * vgas))
            
            if w2 - w1 > 0:
                Z[0] = max(Z)
                if self.CO2_sat:          # CO2 was saturated in previous iteration, but now its not
                    self.CO2_sat = False
                    self.repeat = True
            else:
                Z[0] = min(Z)
                if not self.CO2_sat:
                    self.CO2_sat = True
                    self.repeat = True
        else:
            if self.CO2_sat:          # CO2 was saturated in previous iteration, but now its not
                self.CO2_sat = False
                self.repeat = True
                
            
        return np.real(Z[0])

        
    def MolarVolume(self):
    #=======================================================================
    #  Forms Coefficients of RK-EoS and Solves Cubic for Mixture Molar Volume
    #=======================================================================
    
        RTp = RGASCON * self.tKel / self.pBar
        aT12p = self.aMix / (self.pBar * self.tKel**0.5)
    
        #--Coefficients of the (RK) Cubic Eqn. A-2---------------------------------
        e2 = -RTp                                       # Coeffic of V**2
        e1 = -(RTp * self.bMix - aT12p + self.bMix**2)  # Coeffic of V
        e0 = -aT12p * self.bMix                         # Constant
        
        #--Solve the Cubic------------------------------------------------------
        self.MolarVol = self.cubicSolver(e2, e1, e0)


        
    def fugP(self):
    #=======================================================================
    #  Pressure * Fugacity Coefficient of CO2 or H2O in CO2-Water Mixture
    #=======================================================================
        x, y, kij, aMix, aij = self.x, self.y, self.kij, self.aMix, self.aij
        bMix, b = self.bMix, self.b
        vMol = self.MolarVol           # Always vapor molar volume
        yCO2 = max(min(1,y[0]), 0)
        xCO2 = max(min(1,x[0]), 0)
        self.y = np.array([yCO2, 1.0-yCO2])
        self.x = np.array([xCO2, 1.0-xCO2])
        
        for k in range(2):
            
            t1 = b[k] / bMix * (self.pRT * vMol - 1) 
            t2 = -np.log(self.pRT * max(vMol - bMix,1e-9))
            t3 = sum([y[i] * (aij[i][k] + aij[k][i]) for i in range(2)])
            for i in range(2):
                for j in range(2):
                    t3 -= y[i]**2 * y[j] * (kij[i][j] - kij[j][i]) * (aij[i][i] * aij[j][j])**0.5
            t3 += sum([x[k] * x[i] * (kij[k][i] - kij[i][k]) * (aij[i][i] * aij[k][k])**0.5 for i in range(2)])
            t3 /= aMix
            t3 -= b[k] / bMix
            t4 = (aMix / (bMix * RGASCON * self.tKel ** 1.5)) * np.log(vMol / (vMol + bMix))
         
            #==CO2 Fugacity Coefficient=============================================
            logPhi = t1 + t2 + t3 * t4  # Eq. A-8
      
            #=======================================================================
            #  Pressure * Fugacity Coefficient
            #=======================================================================
            self.fugPi[k] = self.pBar * np.exp(logPhi)
    
    def Ktp(self, K0, Vbar):              # Equation 5
        return K0 * np.exp(self.pRT0 * Vbar)
    
    def K_CO2(self): # Units: bar
    #=======================================================================
    #  CO2 K-value at reservoir Pressure
    #=======================================================================
        if self.low_temp:
            x = [1.189, 1.304e-2, -5.446e-5]
            if self.CO2_sat:
                x = [1.169, 1.368e-2, -5.380e-5] # Liquid CO2 below 31 deg C and above CO2 Psat 
        else:
            x = [1.668, 3.992e-3, -1.156e-5, 1.593e-9]

        K0 = 10**self.FT(self.degC, x)
        
        if self.scaled:
            K0_lt = 10**self.FT(self.degC, [1.189, 1.304e-2, -5.446e-5])
            K0 = self.blended_val(K0_lt, K0)
            
        self.K[0] =  self.Ktp(K0, self.vBar[0])
    
    def K_H2O(self): # Units: bar.mol-1
    #=======================================================================
    #  H2O K-value at reservoir pressure
    #=======================================================================
        if self.low_temp:
            x = [-2.209, 3.097e-2, -1.098e-4, 2.048e-7]
        else:
            x = [-2.1077, 2.8127e-2, -8.4298e-5, 1.4969e-7, -1.1812e-10]
        
        K0 = 10**self.FT(self.degC, x)
        
        if self.scaled:
            K0_lt = 10**self.FT(self.degC, [-2.209, 3.097e-2, -1.098e-4, 2.048e-7])
            K0 = self.blended_val(K0_lt, K0)
            
        self.K[1] =  self.Ktp(K0, self.vBar[1])
    
    def mixMolar(self, x1, P1, P2):
    #=======================================================================
    #  Calculate Mixture Property of two Molar Species
    #=======================================================================
        return x1 * P1 + (1 - x1) * P2
    
    def FT(self, t, x):                        # Equation 6
        return sum([x[i]*t**i for i in range(len(x))])
    
    def calc_gammas(self):
        if self.low_temp:
            return [1,1]
        Am = -3.084e-2*(self.tKel - 373.15) + 1.927e-5*(self.tKel - 373.15)**2                                 # Eq 15
        self.gamma = [np.exp(2*Am*self.x[0]*(1-self.x[0])**2), np.exp((Am - 2*Am*(1-self.x[0]))*self.x[0]**2)] # Eq 12, 13
        
    def A_B(self):
        A = self.K[1] * self.gamma[1] / self.fugPi[1]                                 # Eq 10.
        B = self.fugPi[0] /(CONMOLA * self.gamma[0] * self.gamma_prime * self.K[0])   # Eq 17.
        B = max(EPS,min(1-EPS, B))
        self.A = A
        self.Bprime = B
    
    #============================================================================
    #  A Phase-Partitioning Model for CO2–Brine Mixtures at Elevated Temperatures 
    #  and Pressures: Application to CO2-Enhanced Geothermal Systems
    #  Nicolas Spycher & Karsten Pruess, Transp Porous Med (2010) 82:173–196
    #  DOI 10.1007/s11242-009-9425-y
    #============================================================================
    
    def co2BrineSolubility(self):
        """ Calculates CO2 saturated Brine mutual solubilities
    
            Inputs:
                pBar: Pressure (Bar)
                degC: Temperature (deg C)
                ppm: NaCL equivalent concentration in brine (wt salt per million weight of brine)
    
            Returns tuple of:
                xCO2: Mole fraction CO2 in brine at pBar
                yH2O: Mole fraction H2O in CO2 rich gas
                rhoGas: CO2 rich gas density (gm/cm3)
        """
        pBar = self.pBar
        degC = self.degC
        ppm = self.ppm
        
        #  Initialisation       
        if pBar <= 1.0:
            pBar = 1+EPS
        
        self.tKel = degC + CEL2KEL
        self.pRT = pBar / (RGASCON * self.tKel)

        
        if ppm == None:
            if self.ppm == None: # If salt concentration not defined, set to zero
                self.ppm = 0     # If previously defined, but not redefined here, then leave as it was
        else:
            self.ppm = ppm       # Else if explcitly defined as function parameters, overwrite
        
        if self.degC <= 100:          # Reference pressures delineated by 100 deg C cutoff. 
            self.P0 = 1.0              # Reference Pressure (1 bar at < 100 degC)
        else:                     # Ref Pressure (Water saturation pressure Bar at >= 100 degC)
            self.P0 = self.FT(self.degC, [-1.9906e-1, 2.0471e-3, 1.0152e-4, -1.4234e-6, 1.4168e-8]) 
    
        self.pRT0 = (self.pBar - self.P0) / (RGASCON * self.tKel)
        self.pRT = self.pBar / (RGASCON * self.tKel)
        
        fppM = ppm / 1e6                                                       #--Weight Fraction from PPM
        self.molaL = self.ppm2Molality()                                       #--Molality [gmol/kg]
        self.xSalt = 2 * fppM * MWWAT / (fppM * MWWAT + (1.0 - fppM) * MWSAL)  #--Initial guess of xSalt = xNa + xCl (ie double xNaCl)
        
        self.gammaCO2()  # ----- Calculate Gamma Prime
        
        self.calc_type() # ----- Update the calculation type
        
        #== Component molar volumes
        if self.low_temp:
            VCO2 = 32.6                                 #--CO2 Molar Volume [cm3/gmol]
            VH2O = 18.1                                 #--H2O Molar Volume [cm3/gmol]
        else:
            VCO2 = self.FT(self.tKel - 373.15, [32.6, 3.413e-2])  #--CO2 Molar Volume [cm3/gmol]
            VH2O = self.FT(self.tKel - 373.15, [18.1, 3.137e-2])  #--H2O Molar Volume [cm3/gmol]
        self.vBar = [VCO2, VH2O]
        
        
        #==k-Parameters (y/x partitioning factors) Equations 5 and 6 ==========================
        self.K_H2O()
        i =0
        while not self.repeat:
            i+=1
            if i > 1:
                break
                
            #==k-Parameters (y/x partitioning factors) Equations 5 and 6 ==========================
            self.K_CO2()
            
            #== First estimates of yCO2 and xCO2 ================================================
            if not self.low_temp:
                yCO2 = 1 - self.est_yH2O() # Paper suggests using Psat/P, but this empirical fit appears more accurate
                xCO2 = yCO2 / self.K[0]
            else:
                yCO2 = 1.0
                xCO2 = 0.0
            
            self.y = np.array([yCO2, 1.0-yCO2])
            self.x = np.array([xCO2, 1.0-xCO2])    
            
            # Trigger mixing rules
            self.aMix_RK()
            self.bMix_RK()

            #--Solve the Cubic Equation A-2 --------------------------------------------
            self.MolarVolume()
        
        #--Calculate Fugacity Coefficients*Pressure---------------------------------------
        self.fugP()
        
        if self.scaled: # Recalculate mixing rules and fugacity pressure products with low temp relationships
            phiP_ht = self.fugPi
            self.low_temp = True  # Set flag for low temp coefficiencts
            self.aMix_RK()
            self.bMix_RK()
            
            self.fugP()
            self.fugPi[0] = self.blended_val(self.fugPi[0], phiP_ht[0])
            self.fugPi[1] = self.blended_val(self.fugPi[1], phiP_ht[1])
            self.low_temp = False # Reset low temp flag

        self.calc_gammas()
        self.A_B()     
        
        
        # Calculate results. 
        self.y[1] = (1 - self.Bprime) * CONMOLA / ((1/self.A - self.Bprime) * (2.0 * self.molaL + CONMOLA) + 2 * self.molaL * self.Bprime) # Eq B-7
        
        self.x[0] = self.Bprime * (1.0 - self.y[1])
        self.y[0] = 1.0 - self.y[1]

        mCO2 = self.x[0] * (CONMOLA + 2 * self.molaL) / (1 - self.x[0])      # Eq B-6
        self.xSalt = 2.0 * self.molaL / (2.0 * self.molaL + CONMOLA + mCO2)  # Eq B-3. The 2.0x is stoichiometric ions for NaCl
        self.x[1] = 1.0 - self.x[0] - self.xSalt
        self.x[1] = min(max(self.x[1], 0), 1)

        if not self.low_temp:  # Iterate to solution
            err = 1
            iternum = 0
            while err > EPS and iternum < 100:
                yH2O_last = max(self.y[1], EPS)
                
                # Trigger mixing rules for changed compositions
                self.aMix_RK()
                self.bMix_RK()
        
                #--Fugacity Coefficients*Pressure---------------------------------------
                self.fugP()
                
                if self.scaled:
                    phiP_ht = self.fugPi
                    self.low_temp = True  # Set flag for low temp coefficient calculations
                    self.aMix_RK()
                    self.bMix_RK()
            
                    self.fugP()
                    self.fugPi[0] = self.blended_val(self.fugPi[0], phiP_ht[0])
                    self.fugPi[1] = self.blended_val(self.fugPi[1], phiP_ht[1])
                    self.low_temp = False # Reset low temp flag
                
                self.calc_gammas()
                self.A_B()
                
                self.y[1] = (1 - self.Bprime) * CONMOLA / ((1/self.A - self.Bprime) * (2.0 * self.molaL + CONMOLA) + 2 * self.molaL * self.Bprime) # Eq B-7
                self.y[1] = max(EPS, min(1-EPS, self.y[1]))
                self.x[0] = self.Bprime * (1.0 - self.y[1])
                self.x[0] = max(EPS, min(1-EPS, self.x[0]))
                self.y[0] = 1.0 - self.y[1]
                
                mCO2 = self.x[0] * (CONMOLA + 2 * self.molaL) / (1 - self.x[0])
                
                self.xSalt = 2.0 * self.molaL / (2.0 * self.molaL + CONMOLA + mCO2)  # Eq B-3. The 2.0x is stoichiometric ions for NaCl
                self.x[1] = 1.0 - self.x[0] - self.xSalt
                
                err = abs(self.y[1]/yH2O_last-1)
                
                iternum += 1
        
        #=======================================================================
        #  Re-Compute the CO2/H2O Gas Phase Density
        #=======================================================================
        
        #--Solve the Cubic------------------------------------------------------
        #self.MolarVolume()                                   #--Molar Volume [cm3/gmol]
        self.MwGas = self.mixMolar(self.y[0], MWCO2, MWWAT)   #--Mole weight  [gm/gmol]
        self.rhoGas = self.MwGas / self.MolarVol              #--Density of Gas Mixture [gm/cm3]
        self.GASZ = self.MolarVol * self.pRT
        
    def brine(self):
        results = self.brine_props_co2(self.pBar, self.degC, self.ppm, self.x[0], self.MwBrine)
        self.bDen, self.bVis, self.bVisblty, self.bw, self.Rs, self.Cf_usat = results
        self.Rs -= self.Rs_STD                        # Subtract dissolved CO2 remaining at standard conditions. Note that this is AFTER density and Bw calculations have been conducted
        self.Rs = max(0, self.Rs)
        self.bVisblty = max(1e-12, self.bVisblty)
        
        if not self.metric:
            self.bVisblty = self.bVisblty / BAR2PSI   # CO2 laden viscosibility (1/psi) (at pressures above Psat)
            self.Rs = self.Rs * BBL2CUFT              # Solution GWR CO2 (scf/stb brine)
            self.Cf_usat = self.Cf_usat / BAR2PSI     # Undersaturated brine compressibility (1/psi)
        
    def brine_props_co2(self, pBar, degc, ppm, xCO2, MwB):
        """ Calculates CO2 saturated Brine properties
            1. Pure Brine Density:            Spivey et al. (modified)
            2. Pure Brine viscosity:          Mao-Duan (2009)
            3. CO2 Corrected Brine Density:   Garcia (2001)
            4. CO2 Corrected Brine Viscosity: Islam-Carlson (2012)
            
            pBar: Pressure (Bar)
            degc: Temperature (deg C)
            ppm: Salt weight parts per million brine weight parts
            xCO2: Mole fraction CO2 in brine at pBar
            MwB: MW CO2 free Brine
            
            [sg_CO2_Brine, sg_brine], [cP_CO2_brine, cP_brine], viscosblty, [bw, brine_res_vol], rs, c_usat)
            
            Returns tuple of;
             - Brine Density (gm/cm3) [CO2 Saturated Brine, CO2 Free Brine, Pure Water]
             - viscosity (cP)         [CO2 Saturated Brine, CO2 Free Brine, Pure Water]
             - viscosibility (1/psi))
             - Formation volume factor (res vol/ std vol)  [CO2 Saturated Brine, CO2 Free Brine, Pure Water]
             - Brine compressibility as pressure increases with no change to CO2 saturation (1/Bar)
             - Brine compressibility as pressure decreases with reducing CO2 saturation (1/Bar)
             - Brine solution gas ratio (sm3 CO2 / sm3 brine)
        """
        MwG = MWCO2
    
        wt = ppm / 10000                  # Wt % (0 - 100)
        m = (1000 * (wt / 100) / (MWSAL * (1 - (wt / 100))))  # Molar concentration of NaCl from wt % in gram mol/kg 
        
        Mpa = pBar * 0.1                    # Pressure in mPa
        tKel = degc + CEL2KEL               # Temperature in deg K
        
        def Eq41(t, input_array):  # From McCain Petroleum Reservoir Fluid Properties
            t2 = t / 100
            return (input_array[1] * t2 ** 2 + input_array[2] * t2 + input_array[3]) / (input_array[4] * t2 ** 2 + input_array[5] * t2 + 1)
    
        # Table 4-6 Coefficients
        rhow_t70_arr = [0, -0.127213, 0.645486, 1.03265, -0.070291, 0.639589]
        Ewt_arr = [0, 4.221, -3.478, 6.221, 0.5182, -0.4405]
        Fwt_arr = [0, -11.403, 29.932, 27.952, 0.20684, 0.3768]
        
        # Table 4-7 Coefficients
        Dm2t_arr = [0, -0.00011149, 0.000175105, -0.00043766, 0, 0]
        Dm32t_arr = [0, -0.0008878, -0.0001388, -0.00296318, 0, 0.51103]
        Dm1t_arr = [0, 0.0021466, 0.012427, 0.042648, -0.081009, 0.525417]
        Dm12t_arr = [0, 0.0002356, -0.0003636, -0.0002278, 0, 0]
        
        # Table 4-8 Coefficients
        Emt_arr = [0, 0, 0, 0.1249, 0, 0]
        Fm32t_arr = [0, -0.617, -0.747, -0.4339, 0, 10.26]
        Fm1t_arr = [0, 0, 9.917, 5.1128, 0, 3.892]
        Fm12t_arr = [0, 0.0365, -0.0369, 0, 0, 0]
        
        # Table 4-14 Mao-Duan Coefficients
        d = [0, 2885310, -11072.577, -9.0834095, 0.030925651, -0.0000274071, -1928385.1, 5621.6046, 13.82725, -0.047609523, 0.000035545041]
        
        # Table 4-14 Mao-Duan Coefficients
        a = [-0.21319213, 0.0013651589, -0.0000012191756]
        b = [0.069161945, -0.00027292263, 0.0000002085244]
        c = [-0.0025988855, 0.0000077989227]
        
        # Density of pure water at the reference pressure of 70 MPa, ?w(T, 70 MPa), in g/cm3,
        rhow_t70 = Eq41(degc, rhow_t70_arr)                 # Step 1
        
        Ewt, Fwt = Eq41(degc, Ewt_arr), Eq41(degc, Fwt_arr) # Step 2
        
        Dm2t = Eq41(degc, Dm2t_arr)                         # Step 4
        Dm32t = Eq41(degc, Dm32t_arr)
        Dm1t = Eq41(degc, Dm1t_arr)
        Dm12t = Eq41(degc, Dm12t_arr)
        
        Emt = Eq41(degc, Emt_arr)
        Fm32t = Eq41(degc, Fm32t_arr)
        Fm1t = Eq41(degc, Fm1t_arr)
        Fm12t = Eq41(degc, Fm12t_arr)
        
        # -- CO2-Free Brine Density (gm/cm3)
        def brine_denw(Mpa):
            # cw(T, p), in MPa–1, of pure water at temperature T and pressure p,
            cwtp = (1 / 70) * (1 / (Ewt * (Mpa / 70) + Fwt))        # Eq 4.2
    
            #Density of pure water at temperature T and pressure p.
            Iwt70 = (1 / Ewt) * np.log(abs(Ewt + Fwt))              # Eq 4.3
            Iwtp = (1 / Ewt) * np.log(abs(Ewt * (Mpa / 70) + Fwt))  # Eq 4.4
            rhowtp = rhow_t70 * np.exp(Iwtp - Iwt70)                # Eq 4.5
            
            # Density of brine at temperature T and the reference pressure of 70 MPa
            rhobt70 = (rhow_t70 + Dm2t * m **2 + Dm32t * m ** 1.5 + Dm1t * m + Dm12t * m ** 0.5)  # Eq 4.6
            
            # Brine compressibility coefficients Eb(T,m) and Fb(T,m) from equations (4.7) and (4.8).
            Ebtm = Ewt + Emt * m                                                                  # Eq 4.7
            Fbtm = Fwt + Fm32t * m ** 1.5 + Fm1t * m + Fm12t * m ** 0.5                           # Eq 4.8
            
            # Return methane-free brine density as well as fresh water density in g/cm3
            cbtpm = (1 / 70) * (1 / (Ebtm * (Mpa / 70) + Fbtm))                                   # Eq 4.9 (Step 6)
            Ibt70 = (1 / Ebtm) * np.log(abs(Ebtm + Fbtm))                                         # Eq 4.10 (Step 7)
            Ibtpm = (1 / Ebtm) * np.log(abs(Ebtm * (Mpa / 70) + Fbtm))                            # Eq 4.11
            return rhobt70 * np.exp(Ibtpm - Ibt70), rhowtp                                        # Eq 4.12
        
        # -- CO2 free brine viscosity Mao-Duan (2009) + fresh water viscosity
        def vis_brine(Mpa, rhowtp):
            
            #-- Viscosity of pure water - Eq 4.41 - 4.42
            lnuw_tp = sum([d[i] * np.power(tKel, (i - 3)) for i in range(1, 6)])
            lnuw_tp += sum([rhowtp * (d[i] * np.power(tKel, (i - 8))) for i in range(6, 11)])
            uw_tp = np.exp(lnuw_tp)
    
            #-- Calculate relative viscosity of brine.    Eq 4.43 - 4.47
            AA = a[0] + a[1] * tKel + a[2] * tKel * tKel 
            BB = b[0] + b[1] * tKel + b[2] * tKel * tKel
            CC = c[0] + c[1] * tKel
            lnur_tm = AA * m + BB * m ** 2 + CC * m ** 3
            ur_tm = np.exp(lnur_tm)
            
            # And then brine viscosity in Pa.s, converted to cP
            return (ur_tm * uw_tp * 1000, uw_tp*1000)  # cP - Eq 4.48
        
        def partMolVol(degK):
            #  Partial Molar Volume of dissolved CO2: Garcia Eq (3)
            tC = degK - 273.15
            return 37.51 + tC * (-0.09585 + tC * (0.000874 - tC * 0.0000005044))
        
        # -- Correcting brine density for dissolved CO2, JE Garcia, LBNL Report# 49023, Oct 2011, "Density of Aqueous Solutions of CO2"
        def garciaDensity(rhoBRnoCO2, tKel, pBar, ppm, xCO2, MwB, MwG):
            xNotCO2 = 1.0 - xCO2                         #--Brine (H20+Salt) Mole Fraction
            xRat = xCO2 / xNotCO2                        #--Mole Fraction Ratio, Gas/Brine
            mRat = MwG / MwB                             #--Mole Weight Ratio  , Gas/Brine
            vPhi = partMolVol(tKel)                      #--Apparent Molar Volume of Dissolved CO2
            return (1.0 + mRat * xRat) / (vPhi * xRat / MwB + 1.0 / rhoBRnoCO2) # --Equation 18 of Garcia paper
        
        # Correct CO2 free brine viscosity for dissolved CO2
        # Using approach from "Viscosity Models and Effects of Dissolved CO2", Akand W. Islam and Eric S. Carlson (Jul 2012), Energy Fuels 2012, 26, 8, 5330–5336, https://doi.org/10.1021/ef3006228
        def co2_vis_brine(cP_brine, xCO2):
            # Uses CO2 free brine viscosity (cP) and mole fraction CO2 in brine (xCO2), and returns cP
            return cP_brine * (1 + 4.65 * xCO2**1.0134)
    
        # Density and viscosity at specified pressure & temperature (No CO2)
        sg_brine, rhowtp = brine_denw(Mpa)                        # rhowtp is density of pure water
        cP_brine, cP_freshwater = vis_brine(Mpa, rhowtp)
        
        # Correct for CO2 content
        sg_CO2_Brine = garciaDensity(sg_brine, tKel, pBar, ppm, xCO2, MwB, MwG)
        cP_CO2_brine = co2_vis_brine(cP_brine, xCO2)
        
        # Revaluate at +1 bar for viscosibility calculations
        # Use unchanged xCO2, as viscosibility typically used to characterized UNDERsaturated viscosity behaviour
        sg_brine_, rhowtp_ = brine_denw(Mpa + 0.1)
        cP_brine_, cP_freshwater_ = vis_brine(Mpa + 0.1, rhowtp_)
        sg_CO2_Brine_ = garciaDensity(sg_brine_, tKel, pBar + 1, ppm, xCO2, MwB, MwG)
        cP_CO2_brine_ = co2_vis_brine(cP_brine_, xCO2)
        
        # -- Numerically calculate viscosibility from Mao-Duan (2009) base viscosity with Garcia correction for CO2
        dvdpsi = (cP_CO2_brine_ - cP_CO2_brine)  #-- d[Viscosity/dp [cP/bar]
        viscosblty = dvdpsi * 2 / (cP_CO2_brine + cP_CO2_brine_)  # (1/bar)
        
        # Re-evaluate at reservoir temperature and 1 atmosphere (Rhob_atm)
        Rhob_atm, rhowtp_atm_ = brine_denw(PSTND/10)

        # Re-evaluate at standard conditions
        degc = (60 - 32)/1.8 # 60 deg F
        rhow_t70 = Eq41(degc, rhow_t70_arr)
        Ewt = Eq41(degc, Ewt_arr)
        Fwt = Eq41(degc, Fwt_arr)
        Dm2t = Eq41(degc, Dm2t_arr)
        Dm32t = Eq41(degc, Dm32t_arr)
        Dm1t = Eq41(degc, Dm1t_arr)
        Dm12t = Eq41(degc, Dm12t_arr)
        Emt = Eq41(degc, Emt_arr)
        Fm32t = Eq41(degc, Fm32t_arr)
        Fm1t = Eq41(degc, Fm1t_arr)
        Fm12t = Eq41(degc, Fm12t_arr)
        
        # -- CO2-Free Brine & Freshwater Density at standard conditions (gm/cm3)
        sg_SC_Brine, rhowSC = brine_denw(PSTND/10)
        
        # Calculate mass of 1 sm3 of brine without CO2
        brine_mass = sg_SC_Brine * DENW              # kg brine / sm3 (No CO2)
        brine_moles = brine_mass / MwB               # kg Moles of brine per sm3
        
        # Calculate mass of 1 sm3 of fresh water at reservoir conditions
        water_mass = rhowtp * DENW                   # kg water / sm3 (Freshwater)
        
        # Calculate volume at standard conditions of that much freshwater mass
        sc_volume_freshwater = water_mass / DENW / rhowSC # m3 freshwater
        bw_freshwater = 1/sc_volume_freshwater
        
        # Rearrange: xCO2 = co2_moles / (co2_moles + brine_moles)
        co2_moles = brine_moles * xCO2 / (1 - xCO2)    # kg Moles of dissolved CO2 per sm3
        
        co2_mass = co2_moles * MWCO2                   # kg co2 / sm3 brine
        rs = KGMOL2SM3 * co2_moles                     # sm3 CO2 per sm3 Brine (23.545 m3/kgmol at 60 deg F and 1 atm)
        brine_res_vol = brine_mass / (sg_brine * DENW) # reservoir volume of brine (res m3 No CO2)
        
        # Total mass of brine and CO2, divided by density will yield FVF
        tot_mass = co2_mass + brine_mass
        bw = tot_mass / (sg_CO2_Brine * DENW)
        
        # Undersaturated compressibility = 1/V dV/dP
        c_usat = 1 - sg_CO2_Brine / sg_CO2_Brine_ # 1/Bar    
        
        return ([sg_CO2_Brine, sg_brine, rhowtp], [cP_CO2_brine, cP_brine, cP_freshwater], viscosblty, [bw, brine_res_vol, bw_freshwater], rs, c_usat)