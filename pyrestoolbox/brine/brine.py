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

import numpy as np
import numpy.typing as npt
import pandas as pd

from typing import Tuple
from tabulate import tabulate

import pyrestoolbox.gas as gas # Needed for Z-Factor
from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic
from pyrestoolbox.shared_fns import convert_to_numpy, process_output, halley_solve_cubic
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.constants import R, psc, tsc, degF2R, tscr, scf_per_mol, CUFTperBBL, WDEN, MW_CO2, MW_H2S, MW_N2, MW_AIR, MW_H2
from pyrestoolbox.plyasunov.iapws_if97 import rho_if97 as _rho_if97

def _Eq41(t, input_array):
    """Eq 4.1 from McCain Petroleum Reservoir Fluid Properties"""
    t2 = t / 100
    return (
        input_array[1] * t2 ** 2 + input_array[2] * t2 + input_array[3]
    ) / (input_array[4] * t2 ** 2 + input_array[5] * t2 + 1)

# Spivey coefficient tables (shared by brine_props and brine_props_co2)
_RHOW_T70_ARR = [0, -0.127213, 0.645486, 1.03265, -0.070291, 0.639589]
_EWT_ARR = [0, 4.221, -3.478, 6.221, 0.5182, -0.4405]
_FWT_ARR = [0, -11.403, 29.932, 27.952, 0.20684, 0.3768]
_DM2T_ARR = [0, -0.00011149, 0.000175105, -0.00043766, 0, 0]
_DM32T_ARR = [0, -0.0008878, -0.0001388, -0.00296318, 0, 0.51103]
_DM1T_ARR = [0, 0.0021466, 0.012427, 0.042648, -0.081009, 0.525417]
_DM12T_ARR = [0, 0.0002356, -0.0003636, -0.0002278, 0, 0]
_EMT_ARR = [0, 0, 0, 0.1249, 0, 0]
_FM32T_ARR = [0, -0.617, -0.747, -0.4339, 0, 10.26]
_FM1T_ARR = [0, 0, 9.917, 5.1128, 0, 3.892]
_FM12T_ARR = [0, 0.0365, -0.0369, 0, 0, 0]

def brine_props(p: float, degf: float, wt: float=0, ch4_sat: float=0) -> Tuple:
    """ Calculates Brine properties from modified Spivey Correlation per McCain Petroleum Reservoir Fluid Properties pg 160
        Returns Tuple of (Bw (rb/stb), Density (sg), viscosity (cP), Compressibility (1/psi), Rw GOR (scf/stb))
        p: Pressure (psia)
        degf: Temperature (deg F)
        wt: Salt wt% (0-100)
        ch4_sat: Degree of methane saturation (0 - 1)
    """
    if p <= 0:
        raise ValueError("Pressure must be positive")
    if wt < 0 or wt >= 100:
        raise ValueError(f"Salt weight percent must be >= 0 and < 100, got {wt}")

    Eq41 = _Eq41

    Mpa = p * 0.00689476  # Pressure in mPa
    degc = (degf - 32) / 1.8  # Temperature in deg C
    degk = degc + 273  # Temperature in deg K
    m = (
        1000 * (wt / 100) / (58.4428 * (1 - (wt / 100)))
    )  # Molar concentration of NaCl from wt % in gram mol/kg water

    rhow_t70_arr = _RHOW_T70_ARR
    Ewt_arr = _EWT_ARR
    Fwt_arr = _FWT_ARR
    Dm2t_arr = _DM2T_ARR
    Dm32t_arr = _DM32T_ARR
    Dm1t_arr = _DM1T_ARR
    Dm12t_arr = _DM12T_ARR
    Emt_arr = _EMT_ARR
    Fm32t_arr = _FM32T_ARR
    Fm1t_arr = _FM1T_ARR
    Fm12t_arr = _FM12T_ARR

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
    rhowtp_spivey = rhow_t70 * np.exp(Iwtp - Iwt70)  # Eq 4.5 (Spivey freshwater)
    rhowtp = _rho_if97(degc + 273.15, Mpa) / 1000.0  # IAPWS-IF97 freshwater density (g/cm3)

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
    Rhob_tpm_spivey = rhobt70 * np.exp(
        Ibtpm - Ibt70
    )  # Eq 4.12 - Spivey brine density (no methane)
    # Apply Spivey salt ratio to IAPWS freshwater density
    salt_ratio = Rhob_tpm_spivey / rhowtp_spivey if m > 0 else 1.0
    Rhob_tpm = rhowtp * salt_ratio

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
    rhow_sc_spivey = rhow_sc70 * np.exp(Iw_sc - Iw_sc70)
    rhow_sc = _rho_if97(273.15 + 15, 0.1013) / 1000.0  # IAPWS freshwater at SC
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
    Rhob_scm_spivey = rhob_sc70 * np.exp(
        Ib_scm - Ib_sc70
    )  # Spivey brine density at standard conditions
    salt_ratio_sc = Rhob_scm_spivey / rhow_sc_spivey if m > 0 else 1.0
    Rhob_scm = rhow_sc * salt_ratio_sc

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
    except (ValueError, FloatingPointError):
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

    zee = gas.gas_z(p=p, sg=0.5537, degf=degf, zmethod='BNS',
                    co2=0, h2s=0, n2=0, h2=0)  # Z-Factor of pure methane

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



    zee_sc = gas.gas_z(p=psc, sg=0.5537, degf=tsc, zmethod='BNS',
                       co2=0, h2s=0, n2=0, h2=0)
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
#  A Phase-Partitioning Model for CO2�Brine Mixtures at Elevated Temperatures 
#  and Pressures: Application to CO2-Enhanced Geothermal Systems
#  Nicolas Spycher & Karsten Pruess, Transp Porous Med (2010) 82:173�196
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
                mix = brine.CO2_Brine_Mixture(pres = 5000, temp = 275, ppm = 30000, metric = False)
                mix.bw  # Returns [CO2 Saturated Brine Bw, Pure Brine Bw, Pure Water Bw]
                >> [1.1085795290443725, 1.0543051245909865, 1.0542061001251017]
                
                mix.x  # Returns [xCO2, xBrine]
                >> array([0.02431245, 0.95743175])
                
            Usage example for 175 Bara x 85 degC and 0% NaCl brine:
                mix = brine.CO2_Brine_Mixture(pres = 175, temp = 85)
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

        # Try Halley for all roots
        roots = halley_solve_cubic(e2, e1, e0, flag=0)

        # Fallback to np.roots if Halley returned None
        if roots is None:
            Z = np.roots(np.array([1.0, e2, e1, e0]))
            Z = np.array([x for x in Z if np.isreal(x)])  # Keep only real results
            roots = np.real(Z)

        if len(roots) > 1: # Evaluate which root to use per Eqs 25 and 26 in Spycher & Pruess (2003)
            vgas, vliq = max(roots), min(roots)

            w1 = self.pBar*(vgas - vliq)
            w2 = RGASCON * self.tKel * np.log((vgas - self.bMix)/(vliq - self.bMix)) + self.aMix/(self.tKel**0.5 * self.bMix) * np.log((vgas + self.bMix) * vliq / ((vliq + self.bMix) * vgas))

            if w2 - w1 > 0:
                result = max(roots)
                if self.CO2_sat:          # CO2 was saturated in previous iteration, but now its not
                    self.CO2_sat = False
                    self.repeat = True
            else:
                result = min(roots)
                if not self.CO2_sat:
                    self.CO2_sat = True
                    self.repeat = True
        else:
            result = roots[0]
            if self.CO2_sat:          # CO2 was saturated in previous iteration, but now its not
                self.CO2_sat = False
                self.repeat = True


        return np.real(result)

        
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
    #  A Phase-Partitioning Model for CO2�Brine Mixtures at Elevated Temperatures 
    #  and Pressures: Application to CO2-Enhanced Geothermal Systems
    #  Nicolas Spycher & Karsten Pruess, Transp Porous Med (2010) 82:173�196
    #  DOI 10.1007/s11242-009-9425-y
    #============================================================================
    
    def co2BrineSolubility(self):
        """ Calculates CO2 saturated Brine mutual solubilities
    
            Inputs:
                pBar: Pressure (Bar)
                degC: Temperature (deg C)
                ppm: NaCL equivalent concentration in brine (wt salt per million weight of brine)
    
            Returns Tuple of:
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
            1. Pure Water Density:            IAPWS-IF97 Region 1
            2. Brine Salt Correction:         Spivey et al. (modified)
            3. Pure Brine viscosity:          Mao-Duan (2009)
            4. CO2 Corrected Brine Density:   Garcia (2001)
            5. CO2 Corrected Brine Viscosity: Islam-Carlson (2012)
            
            pBar: Pressure (Bar)
            degc: Temperature (deg C)
            ppm: Salt weight parts per million brine weight parts
            xCO2: Mole fraction CO2 in brine at pBar
            MwB: MW CO2 free Brine
            
            [sg_CO2_Brine, sg_brine], [cP_CO2_brine, cP_brine], viscosblty, [bw, brine_res_vol], rs, c_usat)
            
            Returns Tuple of;
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
        
        Eq41 = _Eq41
        rhow_t70_arr = _RHOW_T70_ARR
        Ewt_arr = _EWT_ARR
        Fwt_arr = _FWT_ARR
        Dm2t_arr = _DM2T_ARR
        Dm32t_arr = _DM32T_ARR
        Dm1t_arr = _DM1T_ARR
        Dm12t_arr = _DM12T_ARR
        Emt_arr = _EMT_ARR
        Fm32t_arr = _FM32T_ARR
        Fm1t_arr = _FM1T_ARR
        Fm12t_arr = _FM12T_ARR
        
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
        def brine_denw(Mpa, tKel_local=tKel):
            # Spivey freshwater density at T, P
            Iwt70 = (1 / Ewt) * np.log(abs(Ewt + Fwt))              # Eq 4.3
            Iwtp = (1 / Ewt) * np.log(abs(Ewt * (Mpa / 70) + Fwt))  # Eq 4.4
            rhowtp_spivey = rhow_t70 * np.exp(Iwtp - Iwt70)         # Eq 4.5
            rhowtp = _rho_if97(tKel_local, Mpa) / 1000.0            # IAPWS freshwater (g/cm3)

            # Spivey brine density at T, P
            rhobt70 = (rhow_t70 + Dm2t * m **2 + Dm32t * m ** 1.5 + Dm1t * m + Dm12t * m ** 0.5)  # Eq 4.6
            Ebtm = Ewt + Emt * m                                                                  # Eq 4.7
            Fbtm = Fwt + Fm32t * m ** 1.5 + Fm1t * m + Fm12t * m ** 0.5                           # Eq 4.8
            Ibt70 = (1 / Ebtm) * np.log(abs(Ebtm + Fbtm))                                         # Eq 4.10
            Ibtpm = (1 / Ebtm) * np.log(abs(Ebtm * (Mpa / 70) + Fbtm))                            # Eq 4.11
            rhob_spivey = rhobt70 * np.exp(Ibtpm - Ibt70)                                         # Eq 4.12

            # Apply Spivey salt ratio to IAPWS freshwater
            salt_ratio = rhob_spivey / rhowtp_spivey if m > 0 else 1.0
            return rhowtp * salt_ratio, rhowtp
        
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
        # Using approach from "Viscosity Models and Effects of Dissolved CO2", Akand W. Islam and Eric S. Carlson (Jul 2012), Energy Fuels 2012, 26, 8, 5330�5336, https://doi.org/10.1021/ef3006228
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
        tKel_sc = (60 - 32) / 1.8 + 273.15  # 60 degF -> K
        sg_SC_Brine, rhowSC = brine_denw(PSTND/10, tKel_local=tKel_sc)
        
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


def make_pvtw_table(
    pi: float,
    degf: float,
    wt: float = 0,
    ch4_sat: float = 0,
    pmin: float = 500,
    pmax: float = 10000,
    nrows: int = 20,
    export: bool = False,
) -> dict:
    """ Generates a PVTW (water PVT) table over a pressure range using brine_props (Spivey correlation).
        Follows the pattern of make_bot_og from the oil module.

        pi: Initial (reference) pressure (psia)
        degf: Temperature (deg F)
        wt: Salt wt% (0-100), default 0
        ch4_sat: Degree of methane saturation (0 - 1), default 0
        pmin: Minimum pressure for table (psia), default 500
        pmax: Maximum pressure for table (psia), default 10000
        nrows: Number of rows in table, default 20
        export: If True, writes PVTW.INC (ECLIPSE keyword) and pvtw_table.xlsx

        Returns dict with keys:
            table: pandas DataFrame with columns Pressure, Bw, Density, Viscosity, Cw, Rsw
            pref: Reference pressure (psia)
            bw_ref: Bw at reference pressure (rb/stb)
            cw_ref: Compressibility at reference pressure (1/psi)
            visw_ref: Viscosity at reference pressure (cP)
            rsw_ref: Rsw at reference pressure (scf/stb)
            den_ref: Density (sg) at reference pressure
    """
    # Build pressure grid, ensuring pi is included
    pressures = list(np.linspace(pmin, pmax, nrows))
    if pi not in pressures:
        pressures.append(pi)
    pressures = sorted(set(pressures))

    bws, dens, visws, cws, rsws = [], [], [], [], []
    for p in pressures:
        bw, lden, visw, cw, rsw = brine_props(p=p, degf=degf, wt=wt, ch4_sat=ch4_sat)
        bws.append(bw)
        dens.append(lden)
        visws.append(visw)
        cws.append(cw)
        rsws.append(rsw)

    df = pd.DataFrame()
    df["Pressure (psia)"] = pressures
    df["Bw (rb/stb)"] = bws
    df["Density (sg)"] = dens
    df["Viscosity (cP)"] = visws
    df["Cw (1/psi)"] = cws
    df["Rsw (scf/stb)"] = rsws

    # Reference properties at pi
    bw_ref, den_ref, visw_ref, cw_ref, rsw_ref = brine_props(
        p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat
    )

    if export:
        # Write ECLIPSE PVTW keyword
        # PVTW format: Pref  Bw  Cw  Visw  Viscosibility
        # Viscosibility set to 0 (constant viscosity assumption)
        pvtw_line = f"  {pi:.1f}  {bw_ref:.6f}  {cw_ref:.6e}  {visw_ref:.4f}  0.0"
        fileout = f"-- Generated by pyResToolbox make_pvtw_table\n"
        fileout += f"-- Temperature: {degf:.1f} deg F, Salt: {wt:.1f} wt%, CH4 sat: {ch4_sat:.2f}\n"
        fileout += f"PVTW\n{pvtw_line} /\n"
        with open("PVTW.INC", "w") as f:
            f.write(fileout)

        # Write full table to Excel
        df.to_excel("pvtw_table.xlsx", index=False)

    return {
        "table": df,
        "pref": pi,
        "bw_ref": bw_ref,
        "cw_ref": cw_ref,
        "visw_ref": visw_ref,
        "rsw_ref": rsw_ref,
        "den_ref": den_ref,
    }


# ============================================================================
# Viscosity correction constants for dissolved gases
# ============================================================================
_IC_A_CO2 = 4.65       # Islam-Carlson (2012) CO2 coefficient
_IC_A_H2S = 1.50       # Calibrated from Murphy & Gaines (1974)
_IC_B = 1.0134         # Islam-Carlson exponent

# Ostermann (1985) SPE 14211 CH4 plateau coefficients
# mu_sat/mu_free = c0 + c1*T + c2*T^2  (T in degF)
# NOTE: SPE 14211 prints c2 as 1.0933e-5 (typo); correct is 1.0933e-6
_OST_C0 = 1.109
_OST_C1 = -5.98e-4
_OST_C2 = 1.0933e-6

# HC component molecular weights for SG-based splitting
_MW_C1 = 16.043
_MW_C2 = 30.069
_MW_C3 = 44.096
_MW_NC4 = 58.122

# Mapping from VLE engine gas keys to plyasunov model gas keys
_VLE_TO_PLYASUNOV = {
    'CH4': 'CH4',
    'C2H6': 'C2H6',
    'C3H8': 'C3H8',
    'nC4H10': 'NC4H10',
    'CO2': 'CO2',
    'N2': 'N2',
    'H2S': 'H2S',
    'H2': 'H2',
}

# Molecular weights for dissolved gas Rs calculations (g/mol)
_SW_GAS_MW = {
    'CH4': 16.0425,
    'C2H6': 30.069,
    'C3H8': 44.096,
    'nC4H10': 58.122,
    'CO2': 44.0095,
    'N2': 28.0134,
    'H2S': 34.081,
    'H2': 2.01588,
}

# Plyasunov model (internal submodule)
from pyrestoolbox.plyasunov import V_phi as _plyasunov_V_phi, gas_mw as _plyasunov_gas_mw

# VLE engine (local copy in brine package)
from pyrestoolbox.brine._lib_vle_engine import calc_gas_brine_equilibrium as _calc_gas_brine_equilibrium


class SoreideWhitson:
    """ Soreide-Whitson (1992) VLE model for multicomponent gas solubility in water/brine,
        with Garcia/Plyasunov density corrections and calibrated viscosity corrections.

        Uses the S&W VLE engine for gas-brine equilibrium calculations, supporting
        multicomponent gas mixtures containing: C1, C2, C3, nC4, CO2, H2S, N2, H2.

        Calculates:
        - Mole fraction of dissolved gas components in aqueous phase
        - Mole fraction of vaporised water in gas phase
        - Water content of gas (stb/mmscf, lb/mmscf)
        - Gas solubility in water (sm3/sm3 or scf/stb) — total and per-gas
        - Gas-saturated brine density via Garcia (2001) Eq. 18 with Plyasunov V_phi
        - Gas-saturated brine viscosity with per-gas corrections (Islam-Carlson, Ostermann, Murphy-Gaines)
        - Brine FVF, compressibility, viscosibility

        Supports fresh and saline water (NaCl equivalent).

        Inputs:
            pres: Pressure (Bar / psia)
            temp: Temperature (deg C / deg F)
            ppm: NaCl equivalent weight concentration in ppm (default 0)
            y_CO2: Mole fraction CO2 in dry gas (default 0)
            y_H2S: Mole fraction H2S in dry gas (default 0)
            y_N2: Mole fraction N2 in dry gas (default 0)
            y_H2: Mole fraction H2 in dry gas (default 0)
            sg: Gas specific gravity — used to estimate HC split among C1-C4 (default 0.65)
            metric: Boolean for units (True=metric, False=oilfield). Default True.
            cw_sat: If True, also calculate saturated compressibility (default False)

        Returns object with following calculated properties:
            .x          : Dict of dissolved gas mole fractions, e.g. {'CO2': 0.024, 'CH4': 0.0015}
            .x_total    : Total dissolved gas mole fraction (sum of all x_i)
            .y          : Dict of gas phase compositions (dry basis, normalized)
            .y_H2O      : Water mole fraction in gas phase
            .water_content : Dict with 'y_H2O', 'stb_mmscf', 'lb_mmscf'
            .bDen       : Brine density (g/cm3) [gas-saturated, gas-free brine, freshwater]
            .bVis       : Brine viscosity (cP) [gas-saturated, gas-free brine, freshwater]
            .bVisblty   : Viscosibility (1/Bar or 1/psi)
            .bw         : Brine FVF (rm3/sm3 or rb/stb) [gas-saturated, gas-free, freshwater]
            .Rs         : Dict of per-gas solution ratios, e.g. {'CO2': sm3/sm3, 'CH4': sm3/sm3}
            .Rs_total   : Total Rs (sum of all per-gas Rs)
            .Cf_usat    : Undersaturated compressibility (1/Bar or 1/psi)
            .Cf_sat     : Saturated compressibility (1/Bar or 1/psi) — only if cw_sat=True
            .MwBrine    : MW of gas-free brine (g/mol)
            .gas_comp   : Normalized gas composition used (including estimated HC split)

        Usage examples:
            # Pure CO2 case, oilfield units
            mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=1.0, metric=False)
            mix.Rs  # Returns per-gas Rs dict, e.g. {'CO2': 15.2}

            # Mixed gas, metric units
            mix = brine.SoreideWhitson(pres=200, temp=80, ppm=10000, y_CO2=0.1, y_H2S=0.05, sg=0.7)
            mix.bDen  # Returns [gas-saturated, gas-free, freshwater] densities

        References:
            Soreide, I. and Whitson, C.H., "Peng-Robinson Predictions for Hydrocarbons,
            CO2, N2, and H2S with Pure Water and NaCl Brine", Fluid Phase Equilibria,
            77, 217-240, 1992.

            Garcia, J.E., "Density of Aqueous Solutions of CO2", LBNL Report 49023, 2001.

            Plyasunov, A.V., Fluid Phase Equilibria (2019, 2020, 2021) — V_phi for
            H2, N2, CH4, CO2, C2H6, C3H8, H2S.

            Islam, A.W. and Carlson, E.S. (2012), Energy & Fuels 26(8), 5330-5336.

            Ostermann, R.D. et al. (1985), SPE 14211.

            Murphy, W.R. and Gaines, T.M. (1974), J. Chem. Eng. Data 19(4), 359-362.
    """

    def __init__(self, pres, temp, ppm=0, y_CO2=0, y_H2S=0, y_N2=0, y_H2=0,
                 sg=0.65, metric=True, cw_sat=False):
        if ppm < 0:
            raise ValueError(f"ppm must be non-negative, got {ppm}")
        if ppm >= 1e6:
            raise ValueError(f"ppm must be less than 1,000,000, got {ppm}")
        if pres <= 0:
            raise ValueError("Pressure must be positive")
        non_hc = y_CO2 + y_H2S + y_N2 + y_H2
        if non_hc > 1.0:
            raise ValueError(f"Sum of non-HC gas fractions ({non_hc}) exceeds 1.0")

        self.metric = metric
        self.ppm = ppm

        # Convert to working units (bar, degC)
        if metric:
            self.pBar = pres
            self.degC = temp
        else:
            self.pBar = pres / BAR2PSI
            self.degC = (temp - 32) / 1.8

        self.degF = self.degC * 1.8 + 32
        self.tKel = self.degC + CEL2KEL
        self.psia = self.pBar * BAR2PSI
        self.Mpa = self.pBar * 0.1
        self.wt_pct = ppm / 10000  # wt% (0-100)
        self.wt_frac = ppm / 1e6   # weight fraction

        # Compute MW of gas-free brine
        xNaCl = (ppm / MWSAL) / ((ppm / MWSAL) + (1e6 - ppm) / MWWAT)
        self.MwBrine = xNaCl * MWSAL + (1 - xNaCl) * MWWAT

        # Estimate HC split from SG and build full gas composition
        self.gas_comp = self._estimate_gas_comp(y_CO2, y_H2S, y_N2, y_H2, sg)

        # Lazy-import VLE engine and Plyasunov model
        self._import_dependencies()

        # Calculate saturated compressibility if requested
        if cw_sat:
            dP_bar = 0.5
            self._calc_properties(self.pBar + dP_bar)
            bw1 = self.bw[0]
            Rs1_total = self.Rs_total
            bDen1 = self.bDen[0]

        # Main calculation at specified pressure
        self._calc_properties(self.pBar)

        if cw_sat:
            bw2 = self.bw[0]
            Rs2_total = self.Rs_total
            # Compute gas Bg at reservoir conditions using gas module Z-factor
            gas_sg = sum(
                self.gas_comp.get(g, 0) * _SW_GAS_MW.get(g, 28.97)
                for g in self.gas_comp
            ) / MW_AIR
            zee = gas.gas_z(p=self.psia, sg=gas_sg, degf=self.degF, zmethod='BNS',
                            co2=self.gas_comp.get('CO2', 0),
                            h2s=self.gas_comp.get('H2S', 0),
                            n2=self.gas_comp.get('N2', 0),
                            h2=self.gas_comp.get('H2', 0))
            Bg = PSTND * zee * self.tKel / (TSTND * self.pBar)  # rm3/sm3
            dBwdP = (bw1 - bw2) / dP_bar
            dRsdP = (Rs1_total - Rs2_total) / dP_bar  # sm3/sm3/bar (internal units)
            self.Cf_sat = (1 / bw2) * (-dBwdP + dRsdP * Bg)
            if not self.metric:
                self.Cf_sat /= BAR2PSI
        else:
            self.Cf_sat = None

        # Unit conversions for oilfield output
        if not self.metric:
            self.bVisblty = self.bVisblty / BAR2PSI
            self.Rs = {k: v * BBL2CUFT for k, v in self.Rs.items()}  # sm3/sm3 -> scf/stb
            self.Rs_total = self.Rs_total * BBL2CUFT
            self.Cf_usat = self.Cf_usat / BAR2PSI

    def _import_dependencies(self):
        """No-op. VLE engine now imported at module level."""
        pass

    @staticmethod
    def _estimate_hc_split(target_mw):
        """Split HC fraction among C1-C4 using constrained exponential decay.

        Uses squared-exponential decay (r²) with bisection to match target MW.
        If unconstrained solution gives C1 < 50%, constrains C1 = 50% and
        distributes remainder among C2-C4 with the same decay pattern.

        Ported from ResToolbox3 _estimateHcSplit().

        Args:
            target_mw: Target molecular weight of HC fraction.

        Returns:
            dict with keys 'CH4', 'C2H6', 'C3H8', 'nC4H10' (mole fractions summing to 1).
        """
        MIN_C1_FRAC = 0.50

        # Clamp to valid range
        mw = max(_MW_C1, min(_MW_NC4, target_mw))

        # Pure methane shortcut
        if mw <= _MW_C1 + 0.1:
            return {'CH4': 1.0, 'C2H6': 0.0, 'C3H8': 0.0, 'nC4H10': 0.0}

        def mw_unconstrained(r):
            r2 = r * r
            r4 = r2 * r2
            r6 = r4 * r2
            s = 1.0 + r2 + r4 + r6
            k = 1.0 / s
            return k * (_MW_C1 + r2 * _MW_C2 + r4 * _MW_C3 + r6 * _MW_NC4)

        def mw_constrained(r):
            r2 = r * r
            r4 = r2 * r2
            heavy_sum = 1.0 + r2 + r4
            heavy_mw = (_MW_C2 + r2 * _MW_C3 + r4 * _MW_NC4) / heavy_sum
            return MIN_C1_FRAC * _MW_C1 + (1.0 - MIN_C1_FRAC) * heavy_mw

        # Try unconstrained solution first — bisect on r in [0, 1]
        r_lo, r_hi = 0.0, 1.0
        for _ in range(50):
            r = (r_lo + r_hi) / 2.0
            mw_calc = mw_unconstrained(r)
            if abs(mw_calc - mw) < 0.001:
                break
            if mw_calc < mw:
                r_lo = r
            else:
                r_hi = r

        # Check if C1 fraction is acceptable
        r2 = r * r
        r4 = r2 * r2
        r6 = r4 * r2
        s = 1.0 + r2 + r4 + r6
        c1_frac = 1.0 / s

        if c1_frac >= MIN_C1_FRAC:
            k = 1.0 / s
            return {
                'CH4': k,
                'C2H6': k * r2,
                'C3H8': k * r4,
                'nC4H10': k * r6,
            }

        # Constrained: fix C1 = MIN_C1_FRAC, distribute C2-C4
        max_mw_constr = MIN_C1_FRAC * _MW_C1 + (1.0 - MIN_C1_FRAC) * _MW_NC4
        mw_target = max(_MW_C1, min(max_mw_constr, mw))

        r_lo, r_hi = 0.0, 10.0
        for _ in range(50):
            r = (r_lo + r_hi) / 2.0
            mw_calc = mw_constrained(r)
            if abs(mw_calc - mw_target) < 0.001:
                break
            if mw_calc < mw_target:
                r_lo = r
            else:
                r_hi = r

        r2 = r * r
        r4 = r2 * r2
        heavy_sum = 1.0 + r2 + r4
        heavy_frac = 1.0 - MIN_C1_FRAC

        return {
            'CH4': MIN_C1_FRAC,
            'C2H6': heavy_frac / heavy_sum,
            'C3H8': heavy_frac * r2 / heavy_sum,
            'nC4H10': heavy_frac * r4 / heavy_sum,
        }

    @staticmethod
    def _estimate_gas_comp(y_CO2, y_H2S, y_N2, y_H2, sg):
        """Estimate full gas composition including HC split from SG.

        The hydrocarbon portion (1 - sum of non-HC) is split among C1-C4
        using constrained exponential decay with r² parameter, matching
        the target HC molecular weight implied by the overall gas SG.

        Returns dict with VLE engine keys (e.g. 'CH4', 'C2H6', 'nC4H10').
        """
        y_hc = 1.0 - y_CO2 - y_H2S - y_N2 - y_H2
        if y_hc < 0:
            raise ValueError(
                f"Non-HC gas fractions sum to {1 - y_hc:.4f}, exceeding 1.0"
            )

        comp = {}
        if y_CO2 > 0:
            comp['CO2'] = y_CO2
        if y_H2S > 0:
            comp['H2S'] = y_H2S
        if y_N2 > 0:
            comp['N2'] = y_N2
        if y_H2 > 0:
            comp['H2'] = y_H2

        if y_hc > 1e-10:
            # Compute apparent HC molecular weight from SG
            mw_gas = sg * MW_AIR
            mw_hc_num = mw_gas - y_CO2 * MW_CO2 - y_H2S * MW_H2S - y_N2 * MW_N2 - y_H2 * MW_H2
            mw_hc = mw_hc_num / y_hc

            # Exponential decay split
            hc_split = SoreideWhitson._estimate_hc_split(mw_hc)
            for gas, frac in hc_split.items():
                if frac > 0:
                    comp[gas] = y_hc * frac

        # Normalize
        total = sum(comp.values())
        if total > 0:
            comp = {k: v / total for k, v in comp.items()}

        return comp

    def _calc_properties(self, pBar):
        """Calculate all brine properties at given pressure (bar)."""
        psia = pBar * BAR2PSI
        degf = self.degF
        Mpa = pBar * 0.1
        tKel = self.tKel
        wt = self.wt_pct  # wt% (0-100)
        salinity_wt_pct = self.ppm / 10000  # Same as wt

        # ================================================================
        # Step 1: VLE — Gas-brine equilibrium
        # ================================================================
        x_gas, water_content = _calc_gas_brine_equilibrium(
            salinity_wt_pct=salinity_wt_pct,
            temperature_F=degf,
            pressure_psia=psia,
            y_CH4=self.gas_comp.get('CH4', 0),
            y_C2H6=self.gas_comp.get('C2H6', 0),
            y_C3H8=self.gas_comp.get('C3H8', 0),
            y_nC4H10=self.gas_comp.get('nC4H10', 0),
            y_CO2=self.gas_comp.get('CO2', 0),
            y_N2=self.gas_comp.get('N2', 0),
            y_H2S=self.gas_comp.get('H2S', 0),
            y_H2=self.gas_comp.get('H2', 0),
            method='flash',
            salinity_method='gamma_phi',
            framework='proposed',
        )

        self.x = x_gas
        self.x_total = sum(x_gas.values())
        self.y = dict(self.gas_comp)  # Normalized dry gas composition
        self.y_H2O = water_content['y_H2O']
        self.water_content = water_content

        # ================================================================
        # Step 2: Base brine properties (gas-free)
        #   Freshwater density: IAPWS-IF97 (international reference standard)
        #   Salt correction: Spivey/McCain ratio (Eq 4.6-4.12 vs Eq 4.1-4.5)
        #   Viscosity: Mao-Duan (2009) via brine_props
        # ================================================================
        # brine_props still needed for viscosity, compressibility, Bw, Rsw
        bw_base, den_base_sg, vis_base_cP, cw_base, rsw_base = brine_props(
            p=psia, degf=degf, wt=wt, ch4_sat=0
        )
        bw_fw, den_fw_sg, vis_fw_cP, cw_fw, rsw_fw = brine_props(
            p=psia, degf=degf, wt=0, ch4_sat=0
        )

        # IAPWS-IF97 freshwater density (kg/m3 -> g/cm3)
        rho_fw_gcc = _rho_if97(tKel, Mpa) / 1000.0
        # Salt correction ratio from Spivey (brine/freshwater)
        salt_ratio = den_base_sg / den_fw_sg if wt > 0 else 1.0
        rho_brine_gcc = rho_fw_gcc * salt_ratio

        # Standard conditions
        p_sc_psia = PSTND * BAR2PSI  # ~14.696 psia
        degf_sc = (TSTND - CEL2KEL) * 1.8 + 32  # ~60 degF
        Mpa_sc = PSTND * 0.1  # ~0.101325 MPa
        tKel_sc = TSTND  # ~288.706 K

        bw_sc, den_sc_sg, _, _, _ = brine_props(
            p=p_sc_psia, degf=degf_sc, wt=wt, ch4_sat=0
        )
        _, den_fw_sc_sg, _, _, _ = brine_props(
            p=p_sc_psia, degf=degf_sc, wt=0, ch4_sat=0
        )

        rho_sc_fw_gcc = _rho_if97(tKel_sc, Mpa_sc) / 1000.0
        salt_ratio_sc = den_sc_sg / den_fw_sc_sg if wt > 0 else 1.0
        rho_sc_brine_gcc = rho_sc_fw_gcc * salt_ratio_sc

        # ================================================================
        # Step 3: Density correction via Garcia Eq. 18 + Plyasunov V_phi
        # ================================================================
        if self.x_total > 0:
            # Compute mole-fraction-weighted effective V_phi and MW
            vphi_eff = 0.0
            mw_eff = 0.0
            for gas_vle, x_i in x_gas.items():
                if x_i <= 0:
                    continue
                yi = x_i / self.x_total  # Fraction among dissolved gases
                gas_ply = _VLE_TO_PLYASUNOV.get(gas_vle, gas_vle.upper())
                vphi_eff += yi * _plyasunov_V_phi(gas_ply, tKel, Mpa)
                mw_eff += yi * _plyasunov_gas_mw(gas_ply)

            # Garcia Eq. 18 in g/cm3:
            # rho = (1 + x2*M2/(M1*x1)) / (x2*V_phi/(M1*x1) + 1/rho1)
            x1 = 1.0 - self.x_total
            M1 = MWWAT
            numerator = 1.0 + self.x_total * mw_eff / (M1 * x1)
            denominator = self.x_total * vphi_eff / (M1 * x1) + 1.0 / rho_brine_gcc
            rho_gas_brine_gcc = numerator / denominator
        else:
            rho_gas_brine_gcc = rho_brine_gcc

        self.bDen = [rho_gas_brine_gcc, rho_brine_gcc, rho_fw_gcc]

        # ================================================================
        # Step 4: Viscosity correction
        # ================================================================
        vis_factor = 1.0
        for gas_vle, x_i in x_gas.items():
            if x_i <= 0:
                continue
            gas_upper = gas_vle.upper()
            if gas_upper == 'CO2':
                vis_factor *= (1.0 + _IC_A_CO2 * x_i ** _IC_B)
            elif gas_upper == 'H2S':
                vis_factor *= (1.0 + _IC_A_H2S * x_i ** _IC_B)
            elif gas_upper == 'CH4':
                ratio = _OST_C0 + _OST_C1 * degf + _OST_C2 * degf ** 2
                vis_factor *= max(ratio, 1.0)
            # C2H6, N2, H2, C3H8, nC4H10: no correction (factor *= 1.0)

        vis_gas_brine = vis_base_cP * vis_factor
        self.bVis = [vis_gas_brine, vis_base_cP, vis_fw_cP]

        # ================================================================
        # Step 5: Viscosibility (numerical derivative at P+1 bar)
        # ================================================================
        bw_p1, den_p1_sg, vis_p1_cP, _, _ = brine_props(
            p=(pBar + 1) * BAR2PSI, degf=degf, wt=wt, ch4_sat=0
        )
        vis_p1_corrected = vis_p1_cP * vis_factor  # Same x_gas (undersaturated)
        dvdp = vis_p1_corrected - vis_gas_brine  # cP/bar
        self.bVisblty = dvdp * 2 / (vis_gas_brine + vis_p1_corrected)  # 1/bar

        # ================================================================
        # Step 6: Bw and Rs
        # ================================================================
        # Mass of 1 sm3 of gas-free brine at standard conditions
        brine_mass = rho_sc_brine_gcc * DENW  # kg/sm3
        brine_moles = brine_mass / self.MwBrine  # kg-mol/sm3

        # Per-gas Rs
        self.Rs = {}
        total_gas_mass = 0.0
        for gas_vle, x_i in x_gas.items():
            if x_i <= 0:
                continue
            gas_moles = brine_moles * x_i / (1.0 - self.x_total)  # kg-mol gas / sm3 brine
            gas_mw_val = _SW_GAS_MW.get(gas_vle, _plyasunov_gas_mw(
                _VLE_TO_PLYASUNOV.get(gas_vle, gas_vle.upper())
            ))
            total_gas_mass += gas_moles * gas_mw_val
            # Rs = moles * molar volume at SC (sm3 gas / sm3 brine)
            self.Rs[gas_vle] = KGMOL2SM3 * gas_moles

        self.Rs_total = sum(self.Rs.values())

        # Bw = total mass / (corrected density * DENW)
        tot_mass = total_gas_mass + brine_mass
        self.bw = [
            tot_mass / (rho_gas_brine_gcc * DENW),   # Gas-saturated
            brine_mass / (rho_brine_gcc * DENW),       # Gas-free brine
            rho_sc_fw_gcc / rho_fw_gcc,                   # Freshwater Bw (sc/res density ratio)
        ]

        # ================================================================
        # Step 7: Undersaturated compressibility
        # ================================================================
        # Recalculate density at P+1 bar with same x_gas (undersaturated)
        Mpa_p1 = (pBar + 1) * 0.1
        # IAPWS freshwater density at P+1
        rho_fw_p1_gcc = _rho_if97(tKel, Mpa_p1) / 1000.0
        # Spivey salt ratio at P+1
        _, den_p1_brine_sg, _, _, _ = brine_props(
            p=(pBar + 1) * BAR2PSI, degf=degf, wt=wt, ch4_sat=0
        )
        _, den_p1_fw_sg, _, _, _ = brine_props(
            p=(pBar + 1) * BAR2PSI, degf=degf, wt=0, ch4_sat=0
        )
        salt_ratio_p1 = den_p1_brine_sg / den_p1_fw_sg if wt > 0 else 1.0
        rho_brine_p1_gcc = rho_fw_p1_gcc * salt_ratio_p1

        if self.x_total > 0:
            # Recompute Garcia density at P+1 with same dissolved gas composition
            vphi_eff_p1 = 0.0
            for gas_vle, x_i in x_gas.items():
                if x_i <= 0:
                    continue
                yi = x_i / self.x_total
                gas_ply = _VLE_TO_PLYASUNOV.get(gas_vle, gas_vle.upper())
                vphi_eff_p1 += yi * _plyasunov_V_phi(gas_ply, tKel, Mpa_p1)

            numerator_p1 = 1.0 + self.x_total * mw_eff / (MWWAT * x1)
            denom_p1 = self.x_total * vphi_eff_p1 / (MWWAT * x1) + 1.0 / rho_brine_p1_gcc
            rho_p1_gcc = numerator_p1 / denom_p1
        else:
            rho_p1_gcc = rho_brine_p1_gcc

        self.Cf_usat = 1.0 - rho_gas_brine_gcc / rho_p1_gcc  # 1/bar