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


# Constants
R = 10.731577089016  # Universal gas constant, ft³·psia/°R·lb.mol
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
MW_AIR = 28.97
MW_H2 = 2.016
