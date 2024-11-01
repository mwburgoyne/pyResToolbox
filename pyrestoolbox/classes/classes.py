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

from enum import Enum

class z_method(Enum):  # Gas Z-Factor calculation model
    DAK = 0
    HY = 1
    WYW = 2
    BUR = 3

class c_method(Enum):  # Gas critical properties calculation method
    PMC = 0
    SUT = 1
    BUR = 2

class pb_method(Enum):  # Bubble point calculation method
    STAN = 0
    VALMC = 1
    VELAR = 2

class rs_method(Enum):  # Oil solution gas calculation method
    VELAR = 0
    STAN = 1
    VALMC = 2

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

