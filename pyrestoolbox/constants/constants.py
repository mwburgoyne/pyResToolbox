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

Physical constants, molecular weights, and unit conversion factors.

Constants
---------
R, psc, tsc, degF2R, tscr, scf_per_mol, CUFTperBBL, WDEN
MW_AIR, MW_CO2, MW_H2S, MW_N2, MW_H2

Unit Conversions (FIELD <-> Eclipse METRIC)
-------------------------------------------
PSI_TO_BAR, BAR_TO_PSI, FT_TO_M, M_TO_FT, IN_TO_MM, MM_TO_IN,
BBL_TO_M3, M3_TO_BBL, MSCF_TO_SM3, SM3_TO_MSCF, MMSCF_TO_SM3, SM3_TO_MMSCF,
STB_TO_SM3, SM3_TO_STB, LBCUFT_TO_KGM3, KGM3_TO_LBCUFT,
SCF_PER_STB_TO_SM3_PER_SM3, SM3_PER_SM3_TO_SCF_PER_STB, and others.

Functions
---------
degf_to_degc    Convert degrees Fahrenheit to Celsius
degc_to_degf    Convert degrees Celsius to Fahrenheit
"""

__all__ = [
    # Physical constants
    'R', 'psc', 'tsc', 'degF2R', 'tscr', 'scf_per_mol', 'CUFTperBBL', 'WDEN',
    'MW_AIR', 'MW_CO2', 'MW_H2S', 'MW_N2', 'MW_H2',
    # Temperature conversion functions
    'degf_to_degc', 'degc_to_degf',
    # Unit conversion constants
    'PSI_TO_BAR', 'BAR_TO_PSI', 'FT_TO_M', 'M_TO_FT', 'IN_TO_MM', 'MM_TO_IN',
    'BBL_TO_M3', 'M3_TO_BBL', 'CUFT_TO_M3', 'M3_TO_CUFT',
    'MSCF_TO_SM3', 'SM3_TO_MSCF', 'MMSCF_TO_SM3', 'SM3_TO_MMSCF',
    'STB_TO_SM3', 'SM3_TO_STB',
    'SCF_PER_STB_TO_SM3_PER_SM3', 'SM3_PER_SM3_TO_SCF_PER_STB',
    'STB_PER_MSCF_TO_SM3_PER_SM3', 'SM3_PER_SM3_TO_STB_PER_MSCF',
    'STB_PER_MMSCF_TO_SM3_PER_SM3', 'SM3_PER_SM3_TO_STB_PER_MMSCF',
    'LBCUFT_TO_KGM3', 'KGM3_TO_LBCUFT',
    'INVPSI_TO_INVBAR', 'INVBAR_TO_INVPSI',
    'PSI2CP_TO_BAR2CP', 'BAR2CP_TO_PSI2CP',
    'PSIFT_TO_BARM', 'BARM_TO_PSIFT',
    'SQFT_TO_SQM', 'SQM_TO_SQFT',
    'D_PER_MSCF_TO_D_PER_SM3', 'D_PER_SM3_TO_D_PER_MSCF',
]


# Constants
R = 10.731577089016  # Universal gas constant, ft��psia/�R�lb.mol
psc = 14.696  # Standard conditions pressure (psia)
tsc = 60  # Standard conditions temperature (deg F)
degF2R = 459.67  # Offset to convert degrees F to degrees Rankine
tscr = tsc + degF2R  # Standard conditions temperature (deg R)
MW_AIR = 28.97  # MW of Air
scf_per_mol = R * tscr / psc  # scf/lb-mol (V = ZnRT/P, Z = 1, n = 1)
CUFTperBBL = 5.61458
WDEN = 62.367 # Water Density lb/cuft

MW_CO2 = 44.01
MW_H2S = 34.082
MW_N2 = 28.014
MW_H2 = 2.016

# ---- Unit conversion constants (FIELD <-> Eclipse METRIC) ----
# Pressure: psia <-> barsa
PSI_TO_BAR = 0.0689475729
BAR_TO_PSI = 1.0 / PSI_TO_BAR

# Temperature: degF <-> degC
def degf_to_degc(degf):
    """Convert degrees Fahrenheit to degrees Celsius."""
    return (degf - 32.0) * 5.0 / 9.0

def degc_to_degf(degc):
    """Convert degrees Celsius to degrees Fahrenheit."""
    return degc * 9.0 / 5.0 + 32.0

# Length: ft <-> m
FT_TO_M = 0.3048
M_TO_FT = 1.0 / FT_TO_M

# Diameter/roughness: inches <-> mm
IN_TO_MM = 25.4
MM_TO_IN = 1.0 / IN_TO_MM

# Volume: bbl <-> m3, cuft <-> m3
BBL_TO_M3 = 0.158987295
M3_TO_BBL = 1.0 / BBL_TO_M3
CUFT_TO_M3 = 0.028316846592
M3_TO_CUFT = 1.0 / CUFT_TO_M3

# Gas rate: Mscf/d <-> sm3/d (1 Mscf = 28.316846592 sm3)
MSCF_TO_SM3 = 28.316846592
SM3_TO_MSCF = 1.0 / MSCF_TO_SM3

# Liquid rate: stb/d <-> sm3/d
STB_TO_SM3 = BBL_TO_M3
SM3_TO_STB = M3_TO_BBL

# GOR: scf/stb <-> sm3/sm3
SCF_PER_STB_TO_SM3_PER_SM3 = CUFT_TO_M3 / BBL_TO_M3  # ~0.17810760667903522
SM3_PER_SM3_TO_SCF_PER_STB = 1.0 / SCF_PER_STB_TO_SM3_PER_SM3

# OGR/CGR: stb/Mscf <-> sm3/sm3 (in Eclipse METRIC, OGR/CGR also sm3/sm3)
STB_PER_MSCF_TO_SM3_PER_SM3 = BBL_TO_M3 / (1000.0 * CUFT_TO_M3)  # ~0.005614583333
SM3_PER_SM3_TO_STB_PER_MSCF = 1.0 / STB_PER_MSCF_TO_SM3_PER_SM3

# Density: lb/cuft <-> kg/m3
LBCUFT_TO_KGM3 = 16.01846337
KGM3_TO_LBCUFT = 1.0 / LBCUFT_TO_KGM3

# Compressibility: 1/psi <-> 1/barsa
INVPSI_TO_INVBAR = BAR_TO_PSI
INVBAR_TO_INVPSI = PSI_TO_BAR

# Pseudopressure: psi^2/cP <-> bar^2/cP
PSI2CP_TO_BAR2CP = PSI_TO_BAR ** 2
BAR2CP_TO_PSI2CP = BAR_TO_PSI ** 2

# Gradient: psi/ft <-> bar/m
PSIFT_TO_BARM = PSI_TO_BAR / FT_TO_M
BARM_TO_PSIFT = 1.0 / PSIFT_TO_BARM

# Area: ft^2 <-> m^2
SQFT_TO_SQM = FT_TO_M ** 2
SQM_TO_SQFT = 1.0 / SQFT_TO_SQM

# Non-Darcy coefficient: day/Mscf <-> day/sm3
D_PER_MSCF_TO_D_PER_SM3 = SM3_TO_MSCF
D_PER_SM3_TO_D_PER_MSCF = MSCF_TO_SM3

# Water content / CGR: stb/MMscf <-> sm3/sm3
STB_PER_MMSCF_TO_SM3_PER_SM3 = BBL_TO_M3 / (1e6 * CUFT_TO_M3)
SM3_PER_SM3_TO_STB_PER_MMSCF = 1.0 / STB_PER_MMSCF_TO_SM3_PER_SM3

# Gas rate: MMscf/d <-> sm3/d (1 MMscf = 28316.846592 sm3)
MMSCF_TO_SM3 = 1000.0 * MSCF_TO_SM3  # 28316.846592
SM3_TO_MMSCF = 1.0 / MMSCF_TO_SM3

# FVF: rb/stb (rcf/scf) <-> rm3/sm3 — both are dimensionless ratios, same numeric value
# (reservoir volume / surface volume in consistent units = same ratio regardless of unit system)
# So no conversion needed for Bo, Bg, Bw expressed as FVF ratios.
# Exception: gas Bg in rcf/scf needs conversion to rm3/sm3 — these are the same ratio, no conversion.
