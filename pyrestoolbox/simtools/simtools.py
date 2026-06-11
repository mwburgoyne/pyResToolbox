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

Backward-compatibility shim for the simtools module.

The implementation now lives in focused sub-files:
    _decks.py       ix_extract_problem_cells, zip_check_sim_deck
    _relperm.py     rel_perm_table, corey, LET, jerauld, fitting helpers
    _aquifer.py     influence_tables, rr_solver
    _vfp.py         make_vfpinj, make_vfpprod and keyword formatters
    _pvt_tables.py  make_bot_og, make_pvtw_table

This module re-exports the full public API so that
`from pyrestoolbox.simtools import simtools` continues to work.
"""

__all__ = [
    'ix_extract_problem_cells', 'rel_perm_table',
    'corey', 'LET', 'jerauld', 'is_let_physical',
    'fit_rel_perm', 'fit_rel_perm_best',
    'influence_tables', 'rr_solver',
    'make_vfpprod', 'make_vfpinj', 'make_bot_og', 'make_pvtw_table',
    'zip_check_sim_deck',
]

# Re-exported for backward compatibility (snapshot registered in
# pyrestoolbox._accelerator.RUST_FLAG_REGISTRY; the live flag used by
# influence_tables is the one in _aquifer.py)
from pyrestoolbox._accelerator import RUST_AVAILABLE, _rust_module

from ._decks import ix_extract_problem_cells, zip_check_sim_deck
from ._relperm import (corey, LET, jerauld, is_let_physical,
                       _normalize_saturation, _apply_kr_model,
                       _build_kr_table, rel_perm_table,
                       fit_rel_perm, fit_rel_perm_best)
from ._aquifer import influence_tables, rr_solver, EPS_T, MAX_ITR
from ._vfp import (_format_vfpinj, _format_vfpprod,
                   make_vfpinj, make_vfpprod)
from ._pvt_tables import make_bot_og, make_pvtw_table
