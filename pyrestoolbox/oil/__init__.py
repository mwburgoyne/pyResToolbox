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
oil_bt              Total two-phase oil FVF (Bo + (Rsi-Rs)*Bg)
oil_sg              Oil specific gravity from API
oil_api             API gravity from specific gravity
oil_ja_sg           Jacoby aromaticity SG
oil_twu_props       Twu critical property correlations
oil_rs_st           Standing Rs correlation
oil_rate_radial     Radial oil flow rate (STB/d)
oil_rate_linear     Linear oil flow rate (STB/d)
oil_harmonize         Harmonize consistent Pb, Rsb, and viscosity
oil_harmonize_pb_rsb  Deprecated wrapper for oil_harmonize (returns 3-tuple)
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
    'oil_co', 'oil_bt', 'oil_sg', 'oil_api', 'oil_ja_sg', 'oil_twu_props', 'oil_rs_st',
    'oil_rate_radial', 'oil_rate_linear', 'oil_harmonize', 'oil_harmonize_pb_rsb',
    'sg_evolved_gas', 'sg_st_gas', 'sgg_wt_avg', 'check_sgs',
    'make_bot_og', 'OilPVT',
    # Enum classes re-exported for convenience (oil.pb_method.STAN, etc.)
    'pb_method', 'rs_method', 'bo_method', 'co_method',
]

# Re-export enum classes
from pyrestoolbox.classes import pb_method, rs_method, bo_method, co_method

# Re-export all public names from sub-modules
from ._utils import oil_sg, oil_api, check_sgs, oil_ja_sg, oil_twu_props
from ._separator import sg_evolved_gas, sg_st_gas, sgg_wt_avg, oil_rs_st
from ._density import oil_deno
from ._correlations import oil_pbub, oil_rs_bub, oil_rs, oil_bo, oil_viso
from ._compressibility import oil_co, oil_bt
from ._rate import oil_rate_radial, oil_rate_linear
from ._harmonize import oil_harmonize, oil_harmonize_pb_rsb
from ._tables import make_bot_og
from ._pvt_class import OilPVT
