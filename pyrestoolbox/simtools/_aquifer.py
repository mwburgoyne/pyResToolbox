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

Aquifer modelling and phase equilibrium: Van Everdingen-Hurst influence
tables (AQUTAB) via inverse Laplace transform, and the Rachford-Rice
flash solver.
"""

from typing import Union, Tuple

import numpy as np
import numpy.typing as npt
from tabulate import tabulate
from ilt import gwr
from mpmath import mp

from pyrestoolbox._accelerator import RUST_AVAILABLE, _rust_module

EPS_T = 1e-15
MAX_ITR = 100


def influence_tables(
    ReDs: list,
    min_td: float = 0.01,
    max_td: float = 200,
    n_incr: int = 20,
    M: int = 7,
    export: bool = False,
) -> Tuple:
    """ Returns a Tuple of;
           1. Dimensionless time list
           2. list of lists of dimensionless pressures at each dimensionless time for each dimensionless radius
        and optionally writes out ECLIPSE styled AQUTAB include file
        Solves Van Everdingin & Hurst Constant Terminal Rate solution via inverse Laplace transform
        Note: This will take ~5 seconds per 20 points to generate

        ReDs: A list of dimensionless radii > 1.0
        min_td: Minimum dimensionless time. Default = 0.01
        max_td: Maximum dimensionless time. Dfeault = 200
        n_incr: Number of increments to split dimensionless time into (log transformed), Default = 20
        M: Laplace inversion accuracy. Higher = more accurate, but more time. Generally 6-12 is good range. Default = 7
        export: Boolean value that controls whether an include file with 'INFLUENCE.INC' name is created. Default: False
    """
    if len(ReDs) == 0:
        raise ValueError("ReDs list cannot be empty")
    if min(ReDs) <=1:
        raise ValueError("ReDs must all be strictly greater than 1.0")

    dtD = np.log(max_td / min_td) / n_incr
    tD = [np.exp(x * dtD + np.log(min_td)) for x in range(n_incr + 1)]
    tD = np.array(tD)

    # Rust fast path: Bessel + GWR entirely in Rust (no Python callbacks)
    if RUST_AVAILABLE:
        pDs = _rust_module.influence_tables_rust(tD.tolist(), list(ReDs), M)
    else:
        # Python fallback via ilt library
        def laplace_Ps(s: float, ReD: float):
            x = mp.sqrt(s)
            # pre-calculate duplicated bessel function for greater efficiency
            i1sReD = mp.besseli(1, x * ReD)
            k1sReD = mp.besselk(1, x * ReD)
            numerator = k1sReD * mp.besseli(0, x) + i1sReD * mp.besselk(0, x)
            denominator = s ** 1.5 * (
                i1sReD * mp.besselk(1, x) - k1sReD * mp.besseli(1, x)
            )
            return numerator / denominator

        pDs = []
        for ReD in ReDs:
            pDs.append(gwr(lambda s: laplace_Ps(s, ReD), tD, M))

    if export:
        inc_out = "---------------------------------------\n"
        inc_out += "AQUTAB\n"
        inc_out += "---------------------------------------\n"
        inc_out += "-- Aquifer Influence Tables...\n"
        inc_out += "--  Based on original work by \n"
        inc_out += "--       Van Everdingin & Hurst\n"
        inc_out += "--\n"
        inc_out += "-- Implemented via Inverse Laplace\n"
        inc_out += "-- transform by pyrestoolbox\n"
        inc_out += "---------------------------------------\n"
        inc_out += "--\n"

        header = ["--", "Table #", "re/rw"]
        table = [["--", 1, "Infinity"]]
        for i, ReD in enumerate(ReDs):
            table.append(["--", i + 2, ReD])
        inc_out += tabulate(table, headers=header) + "\n"
        inc_out += "--\n"
        inc_out += "---------------------------------------\n"
        inc_out += "---- Table #1\n"
        inc_out += "----\n"
        inc_out += "---- Internal Eclipse Table for\n"
        inc_out += "---- Infinite Acting Aquifer\n"
        inc_out += "----\n"

        for i, ReD in enumerate(ReDs):
            inc_out += "---------------------------------------\n"
            inc_out += "----  Table #" + str(i + 2) + "\n"
            inc_out += "----\n"
            inc_out += (
                "----  re/rw = " + str(ReD) + ", No flow outer boundary\n"
            )
            header = ["--", "Table #", "re/rw"]
            table = [["--", "tD", "pD"]]
            for j, td in enumerate(tD):
                table.append(["", td, float(pDs[i][j])])
            inc_out += tabulate(table, headers=header) + "\n"
            inc_out += "/\n"
        inc_out += "---------------------------------------\n"
        inc_out += "-- End of AQUTAB Include File\n"
        inc_out += "---------------------------------------\n"

        with open("INFLUENCE.INC", "w") as text_file:
            text_file.write(inc_out)

    return (tD, pDs)

def rr_solver(
    zi: Union[npt.ArrayLike, list], ki: Union[npt.ArrayLike, list]
) -> Tuple[int, np.ndarray, np.ndarray, float, float]:
    """
    Solves for the root of the Rachford-Rice equation using a method that
    gracefully handles catastrophic roundoff errors.
    The method is outlined in 2022 'Fluid Phase Equilibria' paper by M. Nielsen & H. Lia

    Single-phase feeds (all-liquid or all-vapor) are detected and returned as the
    trivial split (V=0 or V=1) rather than passed to the two-phase solver, which
    is only valid inside the two-phase region.

    Input:
    zi: np.Array of Molar composition (Percent, Fraction or Amounts - will be normalized)
    ki: K-values of the respective molar species.

    Output:
    N_it: Number of iterations required to solve (integer)
    yi: Vapor mole fraction compositions (np.array)
    xi: Liquid mole fraction compositions (np.array)
    V: Vapor molar fraction (float)
    L: Liquid molar fraction (float)
    """
    from pyrestoolbox.brine._lib_vle_engine import rr_solver as _brine_rr_solver

    zi = np.asarray(zi, dtype=float)
    ki = np.asarray(ki, dtype=float)

    if len(zi) == 0 or len(ki) == 0:
        raise ValueError("zi and ki arrays must not be empty")
    if len(zi) != len(ki):
        raise ValueError(f"zi and ki arrays must have the same length (got {len(zi)} and {len(ki)})")

    zi_sum = np.sum(zi)
    if zi_sum <= 0:
        raise ValueError("Sum of zi must be positive")
    if np.any(ki <= 0):
        raise ValueError("All ki values must be positive")

    # The Nielsen-Lia solver is only valid inside the two-phase region; for a
    # single-phase feed it returns physically meaningless inf/nan. Detect the
    # single-phase cases up front and return the trivial split, mirroring the
    # guards in brine._lib_vle_engine.solve_rachford_rice.
    zi_norm = zi / zi_sum
    Km1 = ki - 1.0
    if np.sum(zi_norm * Km1) <= 0:
        # All liquid — V = 0
        xi = zi_norm.copy()
        yi = (ki * zi_norm) / np.sum(ki * zi_norm)
        return 0, yi, xi, 0.0, 1.0
    if np.sum(zi_norm * Km1 / ki) >= 0:
        # All vapor — V = 1
        yi = zi_norm.copy()
        xi = (zi_norm / ki) / np.sum(zi_norm / ki)
        return 0, yi, xi, 1.0, 0.0

    return _brine_rr_solver(zi, ki, tol=EPS_T, max_iter=MAX_ITR)
