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
from typing import Union, List, Tuple

import pandas as pd
from tabulate import tabulate
from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic
from pyrestoolbox.shared_fns import convert_to_numpy, process_output, check_2_inputs, bisect_solve
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.constants import R, psc, tsc, degF2R, tscr, scf_per_mol, CUFTperBBL, WDEN, MW_CO2, MW_H2S, MW_N2, MW_AIR, MW_H2

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
                 'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculating critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
               Defaults to 'PMC'
                S: Skin. Defaults to zero if undefined
        D: Non Darcy Skin Factor (day/mscf). Defaults to zero if undefined
        sg: Gas SG relative to air, Defaults to 0.75 if undefined
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined. If > 0, will change zmethod to 'BUR'
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
    """
    k, h, pr, pwf = (
        np.asarray(k),
        np.asarray(h),
        np.asarray(pr),
        np.asarray(pwf),
    )
    
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR'  
    
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    tc, pc = gas_tc_pc(sg, co2, h2s, n2, h2, cmethod.name, tc, pc)

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
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
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
                 'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
               Defaults to 'PMC' if not specified
                sg: Gas SG relative to air, Defaults to 0.75 if not specified
        degf: Reservoir Temperature (deg F).
        co2: Molar fraction of CO2. Defaults to zero if not specified
        h2s: Molar fraction of H2S. Defaults to zero if not specified
        n2: Molar fraction of Nitrogen. Defaults to zero if not specified
        h2: Molar fraction of Hydrogen. Defaults to zero if not specified
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
    """
    k, area, pr, pwf = (
        np.asarray(k),
        np.asarray(area),
        np.asarray(pr),
        np.asarray(pwf),
    )
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR'  
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])

    tc, pc = gas_tc_pc(sg, co2, h2s, n2, h2, cmethod.name, tc, pc)

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
    tr = degf + degF2R
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
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    cmethod: str = "PMC",
    tc: float = 0,
    pc: float = 0,
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
                 'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
    """
    if tc * pc > 0:  # Critical properties have both been user specified
        return (tc, pc)
    
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
    
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
        sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / MW_AIR) / (
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

    elif (cmethod.name == "BUR"): 
        def _compute_hydrocarbon_critical_properties(hydrocarbon_specific_gravity: float) -> Tuple[float, float]:
            v1, v2, v3 = 0.000229975, 0.186415901, 2.470903632
            offset, pl, vl = -0.032901049, 42.6061669, 3007.108548
            hydrocarbon_molecular_weight = MW_AIR * hydrocarbon_specific_gravity
            vc_on_zc = v1 * hydrocarbon_molecular_weight**2 + v2*hydrocarbon_molecular_weight + v3
            tpc_hc = (offset + vc_on_zc) * vl / (offset + vc_on_zc + pl) 
            ppc_hc = tpc_hc * R / vc_on_zc
            return tpc_hc, ppc_hc
            
        if co2 + h2s + n2 + h2 < 1.0: # If not 100% Inerts, then calculate hydrocarbon MW
            hydrocarbon_specific_gravity = (sg - (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2) / MW_AIR) / (1 - co2 - h2s - n2 - h2)
        else:
            hydrocarbon_specific_gravity = 0.75 # Use default value if 100% inerts to avoid numerical problems
        hydrocarbon_specific_gravity = np.max([0.553779772, hydrocarbon_specific_gravity])  # Methane is lower limit
        
        hydrocarbon_molecular_weight = hydrocarbon_specific_gravity * MW_AIR
        tpc, ppc = (_compute_hydrocarbon_critical_properties(hydrocarbon_specific_gravity))

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
    zmethod: str = "DAK",
    cmethod: str = "PMC",
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
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
                 'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
               'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
               'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
               'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
               Defaults to 'PMC'
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        tc: Critical gas temperature (deg R). Uses cmethod correlation if not specified
        pc: Critical gas pressure (psia). Uses cmethod correlation if not specified
    """
    p, is_list = convert_to_numpy(p)

    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
    elif cmethod == 'BUR':
        zmethod = 'BUR'
    elif zmethod == 'BUR':
        cmethod='BUR'    
        
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
            return z - err/fz_burime(z, rhor, c4), err # Returns updated guess for Z, as well as the error
        
        def fz(z, rhor, c4):                          # The DAK Error function
            return z - (1 + c1*rhor + c2*rhor**2 - c3 * rhor**5 + c4)
    
        def fz_burime(z, rhor, c4):                    # Derivative of the DAK Error function
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
                midz = min(max(0.1, z_bur([pr*pc], tr*tc-degF2R)),3) # First guess using explicit calculation method
                mid_err = 1
                while abs(mid_err) > tol and niter < 100:
                    midz, mid_err = new_dak_z(midz, pr, tr)
                    niter += 1
            else:
                midz = brentq(z_err, bounds[0], bounds[1], args=(pr))  
            return midz
        
        zout = [z_dak_calc(pr, tr) for pr in pprs]
            
        return process_output(zout, is_list)

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

    mws = np.array([44.01, 34.082, 28.014, 2.016, 0])
    tcs = np.array([547.416, 672.120, 227.160, 47.430, 1]) # H2 Tc has been modified
    pcs = np.array([1069.51, 1299.97, 492.84, 187.5300, 1])
    ACF = np.array([0.12256, 0.04916, 0.037, -0.21700, -0.03899])
    VSHIFT = np.array([-0.27593, -0.22896, -0.21066, -0.32400, -0.19076])
    OmegaA = np.array([0.427705, 0.436743, 0.457236, 0.457236, 0.457236])
    OmegaB = np.array([0.0696460, 0.0724373, 0.0777961, 0.0777961, 0.0777961]) 
    VCVIS = np.array([1.46020, 1.46460, 1.35422, 0.67967, 0]) # cuft/lbmol    
        
    # Burgoyne tuned Peng Robinson EOS
    # More information about formulation and applicability can be found here; https://github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor
    def z_bur(psias, degf):
        degR = degf + degF2R

        # Analytic solution for real root(s) of cubic polynomial
        # a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
        # Flag = 1 return Max root, = -1 returns minimum root, = 0 returns all real roots
        def cubic_root(a, flag = 0):
            if a[0] != 1:
                a = np.array(a) / a[0] # Normalize to unity exponent for Z^3
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
                
            if flag == -1:      # Return minimum root
                return min(Zs)
            if flag == 1:       # Return maximum root
                return max(Zs)
            return Zs           # Return all roots
                            
        def calc_bips(hc_mw, degf):                                                                     
            degR = degf + degF2R  
            
            # Hydrocarbon-Inert BIPS (Regressed to Wichert & Synthetic GERG Data)
            # BIP = intcpt + degR_slope/degR + mw_slope * hc_mw
            #                      CO2      H2S        N2        H2 
            intcpts = np.array([0.386557, 0.267007, 0.486589, 0.776917])
            mw_slopes = np.array([-0.00219806, -0.00396541, -0.00316789, 0.0106061])
            degR_slopes = np.array([-158.333, -58.611, -226.239, -474.283])

            hc_bips = list(intcpts + degR_slopes/degR + mw_slopes * hc_mw)
            
            # Inert:Inert BIP Pairs
            #            CO2:H2S       CO2:N2     H2S:N2    CO2:H2  H2S:H2  N2:H2
            inert_bips = [0.0600319, -0.229807, -0.18346, 0.646796, 0.65, 0.369087]
            bips = np.array(hc_bips + inert_bips)
            bip_pairs = [(0, 4), (1, 4), (2, 4), (3, 4), (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]  
            bip_matrix = np.zeros((5, 5))
            
            for i in range(5):
                for j in range(5):
                    for p, pair in enumerate(bip_pairs):
                        if (i, j) == pair:
                            bip_matrix[i, j] = bips[p]
                            bip_matrix[j, i] = bips[p]
                            continue
            return bip_matrix
        
        z = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])
        
        #if tc * pc == 0:  # Critical properties have not been user specified
        #    tc_peng, pc_peng = gas_tc_pc(sg, co2, h2s, n2, h2, cmethod = 'BUR')
        
        tcs[-1], pcs[-1] = tc, pc # Hydrocarbon Tc and Pc from SG using Burgoyne correlation
        trs = degR / tcs
        
        m = 0.37464 + 1.54226 * ACF - 0.26992 * ACF**2
        alpha = (1 + m * (1 - np.sqrt(trs)))**2    
        
        kij = calc_bips(mw_hc, degf)
        
        zout = []
        for psia in psias:
            prs = psia / pcs
            Ai, Bi = OmegaA * alpha * prs / trs**2, OmegaB * prs / trs
            A, B = np.sum(z[:, None] * z * np.sqrt(np.outer(Ai, Ai)) * (1 - kij)), np.sum(z * Bi)
        
            # Coefficients of Cubic: a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
            a = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]
            zout.append(cubic_root(a, flag = 1) - np.sum(z * VSHIFT * Bi)) # Volume translated Z 
        return process_output(zout, is_list) 

    zfuncs = {"DAK": zdak, "HY": z_hy, "WYW": z_wyw, "BUR": z_bur}

    if zmethod.name == 'BUR':
        return zfuncs[zmethod.name](p, degf)
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
    ugz = False
) -> np.ndarray:
    """ Returns Gas Viscosity (cP) or Gas Viscosity * Z-Factor
        Uses Lee, Gonzalez & Eakin (1966) Correlation using equations 2.14-2.17 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        Except if 'BUR' Z-Factor method chosen, in which case tuned LBC viscosity model will be used instead

          p: Gas pressure (psia)
          sg: Gas SG relative to air
          degf: Reservoir Temperature (deg F).
          zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                    defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
          cmethod: Method for calculting critical properties
                   'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                   'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                   Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          zee: Gas Z-Factor. If undefined, will trigger Z-Factor calculation.
          ugz: Boolean flag that if True returns ugZ instead of ug
    """
    p, is_list = convert_to_numpy(p)
    zee, _ = convert_to_numpy(zee)
    
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
    elif cmethod == 'BUR':
        zmethod = 'BUR'
    elif zmethod == 'BUR':
        cmethod='BUR'   
        
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])
    
    t = degf + degF2R
    m = MW_AIR * sg
    
    if not check_2_inputs(zee, p): # Need to calculate Z-Factors if not same length / type as p
         zee = gas_z(p, sg, degf, zmethod, cmethod, co2, h2s, n2, h2, tc, pc)
    
    rho = m * p / (t * zee * R * 62.37)
    
    mws = np.array([44.01, 34.082, 28.014, 2.016, 0])
    tcs = np.array([547.416, 672.120, 227.160, 47.430, 1]) # H2 Tc has been modified
    pcs = np.array([1069.51, 1299.97, 492.84, 187.5300, 1])
    ACF = np.array([0.12256, 0.04916, 0.037, -0.21700, -0.03899])
    VSHIFT = np.array([-0.27593, -0.22896, -0.21066, -0.32400, -0.19076])
    OmegaA = np.array([0.427705, 0.436743, 0.457236, 0.457236, 0.457236])
    OmegaB = np.array([0.0696460, 0.0724373, 0.0777961, 0.0777961, 0.0777961]) 
    VCVIS = np.array([1.46020, 1.46460, 1.35422, 0.67967, 0]) # cuft/lbmol    

    # From https://wiki.whitson.com/bopvt/visc_correlations/
    def lbc(Z, degf, psia, sg, co2=0.0, h2s=0.0, n2=0.0, h2 = 0.0):
        if co2 + h2s + n2 + h2 > 1 or co2 < 0 or h2s < 0 or n2 < 0 or h2 < 0:
            return None
        degR = degf + degF2R
        zi = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])
        if n2 + co2 + h2s + h2 < 1:
            sg_hc = (sg - (co2 * mws[0] + h2s * mws[1] + n2 * mws[2] + h2 * mws[3]) / MW_AIR) / (1 - co2 - h2s - n2 - h2)
        else:
            sg_hc = 0.75 # Irrelevant, since hydrocarbon fraction = 0
        
        sg_hc = max(sg_hc, 0.553779772) # Methane is lower limit
        
        hc_gas_mw = sg_hc * MW_AIR
            
        def vcvis_hc(mw): # Returns hydrocarbon gas VcVis for LBC viscosity calculations   
            return  0.057511062 *  mw + 0.478400158 # ft3/lbmol      
                                                                       
        mws[-1]  = hc_gas_mw       
        tcs[-1], pcs[-1] = gas_tc_pc(hc_gas_mw/MW_AIR, cmethod = 'BUR')
        
        VCVIS[-1] = vcvis_hc(hc_gas_mw)
        degR = degf + degF2R
    
        def stiel_thodos(degR, mws):
            #Calculate the viscosity of a pure component using the Stiel-Thodos correlation.
            Tr = degR / tcs
            ui = []
            Tc = tcs * 5/9 # (deg K)
            Pc = pcs / 14.696
            eta = Tc**(1/6) / (mws**(1/2) * Pc**(2/3)) # Tc and Pc must be in degK and Atm respectively
            
            for i in range(len(Tr)):
                if Tr[i] <= 1.5:
                    ui.append(34e-5 * Tr[i]**0.94 / eta[i])
                else:
                    ui.append(17.78e-5 * (4.58 * Tr[i] - 1.67)**(5/8) / eta[i])
            return np.array(ui)
        
        def u0(zi, ui, mws, Z):  # dilute gas mixture viscosity from Herning and Zippener
            sqrt_mws = np.sqrt(mws)
            return np.sum(zi * ui * sqrt_mws)/np.sum(zi * sqrt_mws)
    
        a = [0.1023, 0.023364, 0.058533, -3.92835e-02,  9.28591e-03] # P3 and P4 have been modified
        # Calculate the viscosity of the mixture using the Lorenz-Bray-Clark method.
        rhoc = 1/np.sum(VCVIS*zi)
        Tc = tcs * 5/9    # (deg K)
        Pc = pcs / 14.696 # (Atm)
    
        eta = np.abs(np.sum(zi*Tc)**(1/6)) / (np.abs(np.sum(zi * mws))**0.5 * np.abs(np.sum(zi * Pc))**(2/3)) # Note 0.5 exponent from original paper
        mw = np.sum(zi * mws)
        rhor = psia / (Z * R * degR) / rhoc
        lhs = a[0] + a[1]*rhor + a[2]*rhor**2 + a[3]*rhor**3 + a[4]*rhor**4
        ui = stiel_thodos(degR, mws)
        vis = (lhs**4 - 1e-4)/eta + u0(zi, ui, mws, Z)
        return process_output(vis, is_list)  
        
    if zmethod.name != 'BUR':
        b = 3.448 + (986.4 / t) + (0.01009 * m)  # 2.16
        c = 2.447 - (0.2224 * b)  # 2.17
        a = ((9.379 + (0.01607 * m)) * np.power(t, 1.5) / (209.2 + (19.26 * m) + t))  # 2.15
        ug = process_output(a * 0.0001 * np.exp(b * np.power(rho, c)), is_list)  # 2.14
    else:
        ug = []
        for psia in p:
            ug.append(lbc(zee, degf, psia, sg, co2, h2s, n2, h2))
        ug = process_output(ug, is_list)
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
) -> np.ndarray:
    """ Returns gas compressibility (1/psi) using the 'DAK' Dranchuk & Abou-Kassem (1975) Z-Factor &
        Critical property correlation values if not explicitly specified
        If h2 > 0, will use the 'BUR' Z-Factor method
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
                   'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                    defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculating Tc and Pc
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                 Defaults to 'PMC
    """
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
        
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
    
    return process_output((vol1 - vol2)/((vol1 + vol2)/2), is_list) # 1/psi

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
) -> np.ndarray:
    """ Returns Bg (gas formation volume factor) for natural gas (rcf/scf)
        p: Gas pressure (psia)
        sg: Gas SG relative to air
        degf: Reservoir Temperature (deg F)
        zmethod: Method for calculating Z-Factor
                 'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'HY' Hall & Yarborough (1973)
                 'WYW' Wang, Ye & Wu (2021)
                 'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                 Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Nitrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified          
    """
    p, is_list = convert_to_numpy(p)
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
        
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
) -> np.ndarray:
    """ Returns gas density for natural gas (lb/cuft)
          p: Gas pressure (psia)
          sg: Gas SG relative to air
          degf: Reservoir Temperature (deg F)
          zmethod: Method for calculating Z-Factor
                   'DAK' Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'HY' Hall & Yarborough (1973)
                   'WYW' Wang, Ye & Wu (2021)
                   'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                   defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
          cmethod: Method for calculting critical properties
                   'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                   'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                   'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                   Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified         
    """
    p, is_list = convert_to_numpy(p)

    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
        
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
        
    return process_output(rhog, is_list)

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
                 'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                 Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          rtol: Relative solution tolerance. Will iterate until abs[(p - poverz * Z)/p] < rtol
    """
    
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
        
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

    return process_output(p, is_list)

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
                 'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                 defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                 Defaults to 'PMC'
          co2: Molar fraction of CO2. Defaults to zero if undefined
          h2s: Molar fraction of H2S. Defaults to zero if undefined
          n2: Molar fraction of Nitrogen. Defaults to zero if undefined
          h2: Molar fraction of Hydrogen. Defaults to zero if undefined
          tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
          pc: Critical gas pressure (psia). Calculates using cmethod if not specified
          rtol: Relative solution tolerance. Will iterate until abs[(grad - calculation)/grad] < rtol
    """

    degR = degf + degF2R

    def grad_err(sg, args):
        grad, p, zmethod, cmethod, tc, pc, co2, h2s, n2, h2 = args
        m = sg * MW_AIR
        zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)
        grad_calc = p * m / (zee * R * degR) / 144
        error = (grad - grad_calc) / grad
        return error

    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
        
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
                   'BUR' Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
                   defaults to 'DAK' if not specified, or to 'BUR' if h2 mole fraction is specified
        cmethod: Method for calculting critical properties
                 'SUT' for Sutton with Wichert & Aziz non-hydrocarbon corrections, or
                 'PMC' for Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
                 'BUR' for Burgoyne method (2024). If h2 > 0, then 'BUR' will be used
                 Defaults to 'PMC'
        co2: Molar fraction of CO2. Defaults to zero if undefined
        h2s: Molar fraction of H2S. Defaults to zero if undefined
        n2: Molar fraction of Nitrogen. Defaults to zero if undefined
        h2: Molar fraction of Hydrogen. Defaults to zero if undefined
        tc: Critical gas temperature (deg R). Calculates using cmethod if not specified
        pc: Critical gas pressure (psia). Calculates using cmethod if not specified
    """
    if h2 > 0:
        cmethod = 'BUR' # The Burgoyne PR EOS method is the only one that can handle Hydrogen
        zmethod = 'BUR' 
        
    zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])
    
    def m_p(p, *args):
        # Pseudo pressure function to be integrated
        degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s, h2 = args
        zee = gas_z(p=p, degf=degf, sg=sg, zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2 = h2, tc=tc, pc=pc)    
        mugz = gas_ug(p, sg, degf, zmethod, cmethod, co2, h2s, n2, h2, tc, pc, zee, ugz=True)  # Gas viscosity z-factor product using a precalculated Z factor
        return 2 * p / (mugz)

    if p1 == p2:
        return 0

    return quad(m_p, p1, p2, args=(degf, sg, zmethod, cmethod, tc, pc, n2, co2, h2s, h2), limit=500)[0]

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
    cond_vol = cgr * CUFTperBBL  # cuft/mmscf surface gas
    cond_sg = oil_sg(api_st)
    cond_mass = cond_vol * (cond_sg * WDEN)  # lb
    surface_gas_moles = 1e6 / scf_per_mol  # lb-moles
    surface_gas_mass = sg_g * MW_AIR * surface_gas_moles  # lb
    cond_mw = 240 - 2.22 * api_st  # lb/lb-moles (Standing correlation)
    cond_moles = cond_mass / cond_mw  # lb-moles
    fws_gas_mass = cond_mass + surface_gas_mass  # lb
    fws_gas_moles = cond_moles + surface_gas_moles  # lb-moles
    fws_gas_mw = fws_gas_mass / fws_gas_moles  # lb/lb-moles
    return fws_gas_mw / MW_AIR

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
                    + (-13064.76 / (t + degF2R))
                    + (-7.3037 * np.log(t + degF2R))
                    + (0.0000012856 * ((t + degF2R) * (t + degF2R)))
                )
            )
            / (p)
            + (np.power(10, ((-3083.87 / (t + degF2R)) + 6.69449)))
        )
        * (1 - (0.00492 * 0) - (0.00017672 * (0 * 0)))
        / 8.32
        / 42
    )
    return content
