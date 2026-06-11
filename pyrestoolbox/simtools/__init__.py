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

Simulation-oriented utilities.

Functions
---------
ix_extract_problem_cells  Extract convergence problem cells from IX PRT files
rel_perm_table            Generate ECLIPSE relative permeability tables (Corey/LET/Jerauld)
corey                     Corey relative permeability model
LET                       LET relative permeability model
jerauld                   Jerauld (Arco) two-parameter relative permeability model
is_let_physical           Check LET curve monotonicity and concavity
fit_rel_perm              Fit Corey/LET/Jerauld model to measured kr data
fit_rel_perm_best         Fit all three models, return best
influence_tables          Van Everdingen-Hurst aquifer influence tables (AQUTAB)
rr_solver                 Rachford-Rice flash calculation solver
make_vfpprod              Generate ECLIPSE VFPPROD lift curve table
make_vfpinj               Generate ECLIPSE VFPINJ lift curve table
make_bot_og               Black oil table generation (PVDO/PVDG/PVTO)
make_pvtw_table           Water PVT table generation (PVTW)
zip_check_sim_deck        Validate and zip simulation deck INCLUDE files
"""

__all__ = [
    'ix_extract_problem_cells', 'rel_perm_table',
    'corey', 'LET', 'jerauld', 'is_let_physical',
    'fit_rel_perm', 'fit_rel_perm_best',
    'influence_tables', 'rr_solver',
    'make_vfpprod', 'make_vfpinj', 'make_bot_og', 'make_pvtw_table',
    'zip_check_sim_deck',
]

# Re-export all public names from sub-modules
from ._decks import ix_extract_problem_cells, zip_check_sim_deck
from ._relperm import (corey, LET, jerauld, is_let_physical,
                       rel_perm_table, fit_rel_perm, fit_rel_perm_best)
from ._aquifer import influence_tables, rr_solver
from ._vfp import make_vfpinj, make_vfpprod
from ._pvt_tables import make_bot_og, make_pvtw_table
