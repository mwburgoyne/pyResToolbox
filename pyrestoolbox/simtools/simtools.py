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

from collections import Counter
import glob
import os
from os.path import exists
import zipfile

from typing import Union, Tuple
import numpy as np
import numpy.typing as npt

import pandas as pd
from tabulate import tabulate
from gwr_inversion import gwr
from mpmath import mp

from pyrestoolbox.classes import (kr_family, kr_table, vlp_method,
                                  z_method, c_method, pb_method, rs_method,
                                  bo_method, uo_method, deno_method, co_method)
from pyrestoolbox.constants import psc
from pyrestoolbox.validate import validate_methods
import pyrestoolbox.gas as gas
import pyrestoolbox.brine as brine
import pyrestoolbox.oil as oil
import pyrestoolbox.oil.oil as _oil_impl

EPS_T = 1e-15
MAX_ITR = 100


def ix_extract_problem_cells(filename: str = "", silent: bool = False) -> list:
    """
    Processes Intersect PRT file to extract convergence issue information
    Prints a summary of worst offenders to terminal (if silent=False), and returns a list
    of sorted dataframes summarising all entities in final convergence row in the PRT files
    List returned is [well_pressure_df, grid_pressure_df, sat_change_df, comp_change_df]
    filename: If empty, will search local directory for PRT file and present list to select from if more than one exists.
              If a filename is furnished, or only one file exists, then no selection will be presented
    silent: False will return only the list of dataframes, with nothing echoed to the terminal
            True will return summary of worst entities to the terminal
    """

    if filename != "":  # A Filename has been provided
        if "PRT" not in filename.upper():
            raise ValueError("File name needs to be an IX print file with .PRT extension: " + filename)

    if filename == "":  # Show selection in local directory
        prt_files = glob.glob("*.PRT", recursive=False)
        if len(prt_files) == 0:
            raise FileNotFoundError("No .PRT files exist in this directory")

        if len(prt_files) > 1:
            table = []
            header = [
                "Index",
                "PRT File Name",
            ]  # Print list of options to select from
            for i in range(len(prt_files)):
                table.append([i, prt_files[i]])
            print(tabulate(table, headers=header))
            print(" ")
            prt_file_idx = int(
                input(
                    "Please choose index of PRT file to parse (0 - "
                    + str(len(prt_files) - 1)
                    + ") :"
                )
            )

            if prt_file_idx not in [i for i in range(0, len(prt_files))]:
                raise ValueError("Index entered outside range permitted")
        else:
            prt_file_idx = 0

        filename = prt_files[prt_file_idx]

    if not silent:
        print("Processing " + filename + "\n")
    file1 = open(filename, "r")  # Kept open for line-by-line reading; closed below
    grab_line1 = False
    grab_line2 = False
    ix_found = False
    max_it = 12
    timesteps = []
    tables = []

    while True:
        line = file1.readline()  # Get next line from file
        # if line is empty, end of file is reached
        if (
            "INTERSECT is a mark of Chevron Corporation, Total S.A. and Schlumberger"
            in line
        ):
            ix_found = True
        if not line:
            break
        if (
            "MaxNewtons                    | Maximum number of nonlinear iterations"
            in line
        ):
            line = line.split("|")
            max_it = int(line[3])
            continue
        if "REPORT   Nonlinear convergence at time" in line:
            table = []
            timesteps.append(line.split()[5])
            grab_line1 = True
            continue
        if grab_line1:
            if "Max" in line:
                grab_line2 = True
                continue
        if grab_line2:
            if "|     |" in line:
                tables.append(table)
                grab_line1, grab_line2 = False, False
                continue
            table.append(line)
    file1.close()

    if not ix_found:
        raise ValueError("Does not appear to be a valid IX PRT file: " + filename)

    # Parse all the last lines in each table
    (
        well_pressures,
        grid_pressures,
        saturations,
        compositions,
        scales,
        balances,
    ) = [[] for x in range(6)]

    for table in tables:
        if len(table) == max_it:
            line = table[-1]
            if "*" not in line:  # If within tolerance, skip
                continue
            line = line.split("|")[2:-1]

            if "*" in line[0]:
                well_pressures.append(line[0].split())
            if "*" in line[1]:
                grid_pressures.append(line[1].split())
            if "*" in line[2]:
                saturations.append(line[2].split())
            if "*" in line[3]:
                compositions.append(line[3].split())
            if "*" in line[4]:
                scales.append(line[4].split())
            if "*" in line[5]:
                balances.append(line[5].split())

    # Summarize bad actors
    def most_frequent(List):
        occurence_count = Counter(List)
        item, count = occurence_count.most_common(1)[0]
        return item, count

    well_pressure_wells = [x[1] for x in well_pressures]
    grid_pressure_locs = [x[1] for x in grid_pressures]
    saturation_locs = [x[1] for x in saturations]
    composition_locs = [x[1] for x in compositions]

    headers = [
        "Issue Type",
        "Total Instances",
        "Most Frequent Actor",
        "Instances",
    ]
    data = [
        well_pressure_wells,
        grid_pressure_locs,
        saturation_locs,
        composition_locs,
    ]
    names = [
        "Well Pressure Change",
        "Grid Pressure Change",
        "Grid Saturation Change",
        "Grid Composition Change",
    ]
    dfs, table, problem_data, problem_data_count = [[] for x in range(4)]
    for d, dat in enumerate(data):
        if len(dat) > 0:
            item, freq = most_frequent(dat)
            problem_data.append(item)
            problem_data_count.append(freq)
        else:
            problem_data.append("None")
            problem_data_count.append(0)
        table.append(
            [names[d], len(dat), problem_data[-1], problem_data_count[-1]]
        )
        dfs.append(pd.DataFrame.from_dict(Counter(dat), orient="index"))

    if not silent:
        print(tabulate(table, headers=headers), "\n")

    for df in dfs:
        try:
            df.columns = ["Count"]
            df.sort_values(by="Count", ascending=False, inplace=True)
        except (ValueError, KeyError):
            pass
    return dfs

def LET(s: np.ndarray, L: float, E: float, T: float) -> np.ndarray:
    """
    Returns LET Relative Permeability curve - Lomeland, F.; Ebeltoft, E.; Thomas, W.H. (2005).

    Input:
    s: Normalized saturation of the phase (np.array)
    L, E, T: Three correlation parameters 'L', 'E', and 'T'.

    Output:
    Relative permeability of the phase (np.array)
    """
    denom = s ** L + E * (1 - s) ** T
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    return np.clip(s ** L / denom, 0.0, 1.0)


def corey(s: np.ndarray, n: float) -> np.ndarray:
    """
    Returns Corey Relative Permeability curve.

    Input:
    s: Normalized saturation of the phase (np.array)
    n: Corey curve exponent.

    Output:
    Relative permeability of the phase (np.array)
    """
    return s ** n


def jerauld(s: np.ndarray, a: float, b: float) -> np.ndarray:
    """
    Returns Jerauld (Arco) Relative Permeability curve.

    Two-parameter model: kr = ((1+b) * S^a) / (1 + b * S^c)
    where c = a * (1 + 1/b).

    Input:
    s: Normalized saturation of the phase (np.array)
    a, b: Jerauld correlation parameters (both must be > 0).

    Output:
    Relative permeability of the phase (np.array)
    """
    b = max(b, 1e-10)
    c = a * (1 + 1 / b)
    denom = 1 + b * s ** c
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    return np.clip(((1 + b) * s ** a) / denom, 0.0, 1.0)


def is_let_physical(s: np.ndarray, L: float, E: float, T: float) -> bool:
    """
    Checks if a LET curve is physical (monotonically increasing with no
    spurious inflection-point reversals in concavity).

    Input:
    s: Normalized saturation array (should be sorted, e.g. np.linspace(0, 1, 100))
    L, E, T: LET correlation parameters.

    Output:
    True if the curve passes physicality checks, False otherwise.
    """
    s = np.sort(s)
    kr = LET(s, L, E, T)
    dkr = np.diff(kr)
    if np.any(dkr < -1e-10):
        return False
    ddkr = np.diff(dkr)
    for i in range(1, len(ddkr)):
        if np.sign(ddkr[i]) - np.sign(ddkr[i - 1]) > 0:
            return False
    return True


def _normalize_saturation(s_actual, s_min, s_max):
    """Normalize actual saturation to [0, 1] range."""
    span = s_max - s_min
    if span <= 0:
        return np.zeros_like(s_actual)
    return np.clip((s_actual - s_min) / span, 0, 1)


def _apply_kr_model(sn, krfamily, n_exp, L, E, T, a_jer, b_jer):
    """Apply the selected kr model to normalized saturation.

    sn: Normalized saturation array [0, 1]
    krfamily: kr_family enum instance
    n_exp: Corey exponent
    L, E, T: LET parameters
    a_jer, b_jer: Jerauld parameters
    Returns: kr array
    """
    fname = krfamily.name
    if fname == "COR":
        return corey(sn, n_exp)
    elif fname == "LET":
        return LET(sn, L, E, T)
    elif fname == "JER":
        return jerauld(sn, a_jer, b_jer)
    else:
        raise ValueError(f"Unknown kr family: {fname}")

def rel_perm_table(
    rows: int,
    krtable: kr_table = kr_table.SWOF,
    krfamily: kr_family = kr_family.COR,
    kromax: float = 1,
    krgmax: float = 1,
    krwmax: float = 1,
    swc: float = 0,
    swcr: float = 0,
    sorg: float = 0,
    sorw: float = 0,
    sgcr: float = 0,
    no: float = 1,
    nw: float = 1,
    ng: float = 1,
    Lw: float = 1,
    Ew: float = 1,
    Tw: float = 1,
    Lo: float = 1,
    Eo: float = 1,
    To: float = 1,
    Lg: float = 1,
    Eg: float = 1,
    Tg: float = 1,
    aw: float = 1,
    bw: float = 1,
    ao: float = 1,
    bo: float = 1,
    ag: float = 1,
    bg: float = 1,
    export: bool = False,
) -> pd.DataFrame:
    """ Returns ECLIPSE styled relative permeability tables
        Users need only define parameters relevant to their table / family selection
        rows: Integer value specifying the number of table rows desired
        krtable: A string or kr_table Enum class that specifies one of three table type choices;
                   SWOF: Water / Oil table
                   SGOF: Gas / Oil table
                   SGFN: Gas / Water table
        krfamily: A string or kr_family Enum class that specifies one of three curve function choices;
                   COR: Corey Curve function
                   LET: LET Relative permeability function
                   JER: Jerauld (Arco) two-parameter model
        kromax: Max Kr relative to oil. Default value = 1
        krgmax: Max Kr relative to gas. Default value = 1
        krwmax: Max Kr relative to water. Default value = 1
        swc: Minimum water saturation. Default value = 0
        swcr: Maximum water saturation for imobile water. Default value = 0
        sorg: Maximum oil saturation relative to gas for imobile oil. Default value = 0
        sorw: Maximum oil saturation relative to water for imobile oil. Default value = 0
        sgcr: Maximum gas saturation relative to water for imobile gas. Default value = 0
        no, nw, ng: Corey exponents to oil, water and gas respectively. Default values = 1
        Lw, Ew, Tw: LET exponents to water. Default values = 1
        Lo, Eo, To: LET exponents to oil. Default values = 1
        Lg, Eg, Tg: LET exponents to gas. Default values = 1
        aw, bw: Jerauld parameters for water. Default values = 1
        ao, bo: Jerauld parameters for oil. Default values = 1
        ag, bg: Jerauld parameters for gas. Default values = 1
        export: Boolean value that controls whether an include file with same name as krtable is created. Default: False
    """

    if isinstance(krtable, str):
        try:
            krtable = kr_table[krtable.upper()]
        except KeyError:
            raise ValueError(f"Incorrect table type specified: '{krtable}'")
    if isinstance(krfamily, str):
        try:
            krfamily = kr_family[krfamily.upper()]
        except KeyError:
            raise ValueError(f"Incorrect krfamily specified: '{krfamily}'")

    if rows < 2:
        raise ValueError(f"rows must be at least 2, got {rows}")

    # Consistency checks
    errors = []
    swcr = max(swc, swcr)
    if sorg + sgcr + swc >= 1:
        errors.append("sorg+sgcr+swc must be less than 1")
    if sorw + swcr >= 1:
        errors.append("sorw+swcr must be less than 1")
    if errors:
        raise ValueError("Saturation consistency check failure: " + "; ".join(errors))

    fname = krfamily.name

    # Validate that required parameters are available for the chosen family
    if fname == "COR":
        if krtable.name == "SWOF" and not (no > 0 and nw > 0):
            raise ValueError("Corey SWOF requires no > 0 and nw > 0")
        if krtable.name == "SGOF" and not (no > 0 and ng > 0):
            raise ValueError("Corey SGOF requires no > 0 and ng > 0")
        if krtable.name == "SGWFN" and not (ng > 0 and nw > 0):
            raise ValueError("Corey SGWFN requires ng > 0 and nw > 0")
    elif fname == "LET":
        if krtable.name == "SWOF" and not (Lw > 0 and Ew > 0 and Tw > 0 and Lo > 0 and Eo > 0 and To > 0):
            raise ValueError("LET SWOF requires Lw, Ew, Tw, Lo, Eo, To all > 0")
        if krtable.name == "SGOF" and not (Lg > 0 and Eg > 0 and Tg > 0 and Lo > 0 and Eo > 0 and To > 0):
            raise ValueError("LET SGOF requires Lg, Eg, Tg, Lo, Eo, To all > 0")
        if krtable.name == "SGWFN" and not (Lw > 0 and Ew > 0 and Tw > 0 and Lg > 0 and Eg > 0 and Tg > 0):
            raise ValueError("LET SGWFN requires Lw, Ew, Tw, Lg, Eg, Tg all > 0")
    elif fname == "JER":
        if krtable.name == "SWOF" and not (ao > 0 and bo > 0 and aw > 0 and bw > 0):
            raise ValueError("Jerauld SWOF requires ao, bo, aw, bw all > 0")
        if krtable.name == "SGOF" and not (ao > 0 and bo > 0 and ag > 0 and bg > 0):
            raise ValueError("Jerauld SGOF requires ao, bo, ag, bg all > 0")
        if krtable.name == "SGWFN" and not (ag > 0 and bg > 0 and aw > 0 and bw > 0):
            raise ValueError("Jerauld SGWFN requires ag, bg, aw, bw all > 0")

    if krtable.name == "SWOF":
        # Build saturation array, adjusting ndiv to account for special endpoints
        ndiv = rows
        if swcr > swc:
            ndiv -= 2
        if sorw > 0:
            ndiv -= 1
        ndiv = min(ndiv, rows - 1)
        sw_eps = [swc, swcr, 1 - sorw, 1]
        swn = np.arange(0, 1, 1 / ndiv)
        sw = swn * (1 - swcr - sorw) + swcr
        sw = sorted(set(list(sw) + sw_eps))
        sw = np.array(sw)

        # Water kr: normalize over [swcr, 1-sorw]
        swn_w = _normalize_saturation(sw, swcr, 1 - sorw)
        krw = krwmax * _apply_kr_model(swn_w, krfamily, nw, Lw, Ew, Tw, aw, bw)

        # Oil kr: normalize over [swc, 1-sorw]
        swn_o = _normalize_saturation(sw, swc, 1 - sorw)
        kro = kromax * _apply_kr_model(1 - swn_o, krfamily, no, Lo, Eo, To, ao, bo)

        kr_df = pd.DataFrame()
        kr_df["Sw"] = sw
        kr_df["Krwo"] = krw
        kr_df["Krow"] = kro
        if export:
            df = kr_df.set_index("Sw")
            headings = ["-- Sw", "Krwo", "Krow"]
            fileout = "SWOF\n" + tabulate(df, headings) + "\n/"
            with open("SWOF.INC", "w") as text_file:
                text_file.write(fileout)
        return kr_df

    elif krtable.name == "SGOF":
        # Build saturation array
        ndiv = rows
        if sorg > 0:
            ndiv -= 1
        ndiv = min(ndiv, rows - 1)
        sg_eps = [0, 1 - swc - sorg]
        sgn = np.arange(0, 1, 1 / ndiv)
        sg = sgn * (1 - swc - sorg)
        sg = sorted(set(list(sg) + sg_eps))
        sg = np.array(sg)

        # Gas kr: normalize over [0, 1-swc-sorg]
        sgn_g = _normalize_saturation(sg, 0, 1 - swc - sorg)
        krg = krgmax * _apply_kr_model(sgn_g, krfamily, ng, Lg, Eg, Tg, ag, bg)

        # Oil kr: same normalization, reversed
        kro = kromax * _apply_kr_model(1 - sgn_g, krfamily, no, Lo, Eo, To, ao, bo)

        kr_df = pd.DataFrame()
        kr_df["Sg"] = sg
        kr_df["Krgo"] = krg
        kr_df["Krog"] = kro
        if export:
            headings = ["-- Sg", "Krgo", "Krog"]
            df = kr_df.set_index("Sg")
            fileout = "SGOF\n" + tabulate(df, headings) + "\n/"
            with open("SGOF.INC", "w") as text_file:
                text_file.write(fileout)
        return kr_df

    elif krtable.name == "SGWFN":
        # Build saturation array
        ndiv = rows
        if sgcr > 0:
            ndiv -= 1
        ndiv = min(ndiv, rows - 1)
        sg_eps = [0, sgcr, 1 - swc]
        sgn = np.arange(0, 1, 1 / ndiv)
        sg = sgn * (1 - swc) + 0  # Range is [0, 1-swc]
        sg = sorted(set(list(sg) + sg_eps))
        sg = np.array(sg)

        # Gas kr: normalize over [sgcr, 1-swc]
        sgn_g = _normalize_saturation(sg, sgcr, 1 - swc)
        krg = krgmax * _apply_kr_model(sgn_g, krfamily, ng, Lg, Eg, Tg, ag, bg)

        # Water kr: normalize over [0, 1-swc]
        sgn_w = _normalize_saturation(sg, 0, 1 - swc)
        krw = krwmax * _apply_kr_model(1 - sgn_w, krfamily, nw, Lw, Ew, Tw, aw, bw)

        kr_df = pd.DataFrame()
        kr_df["Sg"] = sg
        kr_df["Krgw"] = krg
        kr_df["Krwg"] = krw
        if export:
            df = kr_df.set_index("Sg")
            headings = ["-- Sg", "Krgw", "Krwg"]
            fileout = "SGWFN\n" + tabulate(df, headings) + "\n/"
            with open("SGWFN.INC", "w") as text_file:
                text_file.write(fileout)
        return kr_df

    raise ValueError("Check that you have specified table type as SWOF, SGOF or SGWFN")


def fit_rel_perm(
    sw: Union[npt.ArrayLike, list],
    kr: Union[npt.ArrayLike, list],
    krfamily: kr_family = kr_family.COR,
    krmax: float = 1.0,
    sw_min: float = 0.0,
    sw_max: float = 1.0,
) -> dict:
    """ Fits a relative permeability model (Corey, LET, or Jerauld) to measured kr data.

        Uses scipy.optimize.least_squares to find the best-fit parameters by minimizing
        the sum of squared residuals between the model and measured data.

        sw: Array of water saturation (or gas saturation) values corresponding to kr data.
        kr: Array of measured relative permeability values.
        krfamily: A string or kr_family Enum class specifying the model to fit;
                  COR: Corey (1 parameter: n)
                  LET: LET (3 parameters: L, E, T)
                  JER: Jerauld (2 parameters: a, b)
        krmax: Maximum kr endpoint value. The model is scaled by this value. Default = 1.0
        sw_min: Minimum saturation endpoint (critical saturation). Used to normalize sw to [0,1]. Default = 0.0
        sw_max: Maximum saturation endpoint. Used to normalize sw to [0,1]. Default = 1.0

        Returns: Dictionary with keys:
            'params': Dict of fitted parameter names and values (e.g. {'n': 2.5} for Corey)
            'krmax': The krmax used
            'sw_min': The sw_min used
            'sw_max': The sw_max used
            'residuals': Array of residuals at the solution
            'ssq': Sum of squared residuals
            'kr_pred': Array of predicted kr values at the input sw points
    """
    from scipy.optimize import least_squares

    sw = np.asarray(sw, dtype=float)
    kr = np.asarray(kr, dtype=float)

    if len(sw) != len(kr):
        raise ValueError(f"sw and kr must have same length (got {len(sw)} and {len(kr)})")
    if len(sw) < 2:
        raise ValueError("Need at least 2 data points for fitting")

    if isinstance(krfamily, str):
        try:
            krfamily = kr_family[krfamily.upper()]
        except KeyError:
            raise ValueError(f"Unknown kr family: '{krfamily}'")

    # Normalize saturation
    sn = _normalize_saturation(sw, sw_min, sw_max)

    # Scale measured kr by krmax
    kr_scaled = kr / max(krmax, 1e-30)

    fname = krfamily.name

    if fname == "COR":
        def residuals(p):
            n = p[0]
            return corey(sn, n) - kr_scaled
        x0 = [2.0]
        bounds = ([0.1], [20.0])
    elif fname == "LET":
        def residuals(p):
            L, E, T = p
            return LET(sn, L, E, T) - kr_scaled
        x0 = [2.0, 1.0, 2.0]
        bounds = ([0.01, 0.01, 0.01], [50.0, 50.0, 50.0])
    elif fname == "JER":
        def residuals(p):
            a, b = p
            return jerauld(sn, a, b) - kr_scaled
        x0 = [2.0, 1.0]
        bounds = ([0.1, 0.01], [20.0, 100.0])
    else:
        raise ValueError(f"Unknown kr family: {fname}")

    result = least_squares(residuals, x0, bounds=bounds)

    # Build prediction at input points
    kr_pred = krmax * (result.fun + kr_scaled)

    # Package results
    if fname == "COR":
        params = {"n": result.x[0]}
    elif fname == "LET":
        params = {"L": result.x[0], "E": result.x[1], "T": result.x[2]}
    elif fname == "JER":
        params = {"a": result.x[0], "b": result.x[1]}

    return {
        "params": params,
        "krmax": krmax,
        "sw_min": sw_min,
        "sw_max": sw_max,
        "residuals": result.fun,
        "ssq": float(np.sum(result.fun ** 2)),
        "kr_pred": kr_pred,
    }


def fit_rel_perm_best(
    sw: Union[npt.ArrayLike, list],
    kr: Union[npt.ArrayLike, list],
    krmax: float = 1.0,
    sw_min: float = 0.0,
    sw_max: float = 1.0,
) -> dict:
    """ Fits all three relative permeability models (Corey, LET, Jerauld) to measured data
        and returns the best fit (lowest sum of squared residuals).

        sw: Array of water saturation (or gas saturation) values corresponding to kr data.
        kr: Array of measured relative permeability values.
        krmax: Maximum kr endpoint value. Default = 1.0
        sw_min: Minimum saturation endpoint (critical saturation). Default = 0.0
        sw_max: Maximum saturation endpoint. Default = 1.0

        Returns: Dictionary with same keys as fit_rel_perm, plus:
            'family': Name of the best-fit model ('COR', 'LET', or 'JER')
            'all_results': Dict of all three fit results keyed by family name
    """
    all_results = {}
    for fam in [kr_family.COR, kr_family.LET, kr_family.JER]:
        all_results[fam.name] = fit_rel_perm(sw, kr, krfamily=fam,
                                              krmax=krmax, sw_min=sw_min, sw_max=sw_max)

    best_name = min(all_results, key=lambda k: all_results[k]["ssq"])
    best = all_results[best_name].copy()
    best["family"] = best_name
    best["all_results"] = all_results
    return best


def influence_tables(
    ReDs: list,
    min_td: float = 0.01,
    max_td: float = 200,
    n_incr: int = 20,
    M: int = 8,
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
        M: Laplace invesrion accuracy. Higher = more accurate, but more time. Generally 6-12 is good range. Default = 8
        export: Boolean value that controls whether an include file with 'INFLUENCE.INC' name is created. Default: False
    """
    if len(ReDs) == 0:
        raise ValueError("ReDs list cannot be empty")
    if min(ReDs) <=1:
        raise ValueError("ReDs must all be strictly greater than 1.0")

    # Eq 23 from SPE 81428
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

    dtD = np.log(max_td / min_td) / n_incr
    tD = [np.exp(x * dtD + np.log(min_td)) for x in range(n_incr + 1)]
    tD = np.array(tD)

    pDs = []
    for ReD in ReDs:
        print("Calculating ReD = " + str(ReD))
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

def zip_check_sim_deck(files2scrape = [], tozip = True, console_summary = True):
    """ Performs recursive ECL/IX deck zip/check
        Crawls through all INCLUDE files in a deck, including an unlimited number of subdirectories and nested INCLUDE references, 
        and (a) checks that all include files exist, then optionally (b) creates a zip file of all associated files
        It does NOT zip any files that are in a higher directory than the .DATA file, but it does flag any such files so users can manually include them
        
        Run in directory containing the simulation DATA or AFI files (or change directory to there)
        Select one or more decks by their index number when prompted, then follow instructions
        
        files2scrape: A list of file names to scrape. These must be .DATA and/or .AFI files
                      If no (or empty) list is passed, then user will be prompted to select one or more from the files existing in the current directory
        tozip: Controls whether all files will be zipped, or whether a check that all files exist is all that is performed.
               If no (or empty) files2scrape list was passed, then user will be prompted for their choice no matter what was specified.   
               
        console_summary: Controls whether summary results are returned to console or not
                         If False, will return a list of missing INCLUDE files. A zero length list would indicate no missing files.
    """

    def get_loc(mask):
        input_files = []
        
        for files in mask:
    
            for file in glob.glob(files.upper()): # Grab a list of all the files in the directory with file mask
                input_files.append(file)
            for file in glob.glob(files.lower()): # In case lowercase extension used
                input_files.append(file)
        input_files = list(set(input_files)) # In case duplicated file names due to checking upper and lower case
        
        if len(input_files) == 0:
            raise FileNotFoundError('No '+mask+' files exist in this directory')
    
        print(' ')
        
        input_files.sort()
        
        table = []
        header=['Index', 'File Name']  # Print list of options to select from
        for i in range(len(input_files)):
            table.append([i,input_files[i]])
        print(tabulate(table,headers=header))    
        print(' ')
        file_idx = input('Please choose index(s) of file to parse separated by commas (0 - '+str(len(input_files)-1)+') :')
        file_idxs = [int(x) for x in file_idx.split(',')]
    
        if not all(item in [i for i in range(0, len(input_files))] for item in file_idxs):
            raise ValueError('Index entered outside range permitted')
        
        in_files =  [input_files[x] for x in file_idxs]
        return in_files
    
    types = ('*.DATA', '*.afi')
    if len(files2scrape) == 0:
        files2scrape = get_loc(types)
        tozip = True
        method = input('Zip or Check files? (Z/c): ')
        if method.upper()=='C':
            tozip = False
    
    files2scrape = [x for x in files2scrape]
    parent_filenames = [x for x in files2scrape]
    files2scrape_set = set(files2scrape)
    missing_parents = []
    # Start stepping through, looking for INCLUDE file statements
    get_include = False
    got_all = False

    nscraped = 0
    higher_dir = False
    missing = []
    while not got_all:
        if console_summary:
            print('Scanning through: '+files2scrape[nscraped])

        # Load file into list
        try:
            lines = list(open(files2scrape[nscraped], 'r'))
        except (FileNotFoundError, PermissionError, OSError):
            if not exists(files2scrape[nscraped]):
                if files2scrape[nscraped] not in missing:
                    missing.append(files2scrape[nscraped])
                    missing_parents.append(parent_filenames[nscraped])
            nscraped += 1
            got_next = True
            if len(files2scrape) == nscraped:
                got_all = True
            continue

        for line in lines:
            line = line.strip() # Remove leading and trailing spaces
            if line[:3]== '--' or line[:1]=='#': # Skip all comments
                continue
            if line[:4]== 'END': # Skip everything after an END command
                break
            if line.upper()[:7]=='INCLUDE':
                get_include = True

            if get_include:

                # Remove any trailing comments that might give false negatives
                line = line.split('--')[0]
                line = line.split('#')[0]

                line = line.replace('"',"'") # Remove double inverted commas, replace with single commas
                if '.' not in line and "'" not in line: # Filename not in this line
                    continue

                parts = line.split("'")
                if len(parts) < 2: # Malformed INCLUDE line
                    get_include = False
                    continue
                include_file = parts[1].strip()
                if include_file not in files2scrape_set:
                    files2scrape.append(include_file)
                    files2scrape_set.add(include_file)
                    parent_filenames.append(files2scrape[nscraped])
                    if '..' in include_file:
                        higher_dir = True
                get_include = False
                continue
 
        got_next = True
                
        # Finished scraping a file - are there any left?
        nscraped += 1
        if len(files2scrape) == nscraped:
            got_all = True
            continue
            
        
        
    if len(missing)>0:
        if console_summary:
            print('\n****** MISSING FILES ******\n')
            for f, file in enumerate(missing):
                print(file, 'from', missing_parents[f])
    
    higher_files = []
    if not tozip:
        if len(missing) == 0:
            if console_summary:
                print('\nALL INCLUDE FILES FOUND\n')
            for file in files2scrape:
                if '..' in file:
                    higher_files.append(file)
            if len(higher_files) > 0:
                if console_summary:
                    print('\n********** WARNING: Some INCLUDE files in a parent directory **********')
                    for file in higher_files:
                            print(file)
                    print('\nThese would need to be added manually to a zip file')
    
    if tozip:
        if console_summary:
            print('\n'+str(len(files2scrape))+' files to zip')
        if len(missing)>0:
            if console_summary:
                cont = input('Continue to zip even with missing files? (Y/n): ')
                if cont.upper() =='N':
                    raise RuntimeError("Zip aborted by user due to missing files")
        lista_files = files2scrape
        dotindex = ''.join(lista_files[0]).rindex('.')
        zipname = lista_files[0][:dotindex]+'.zip'
        with zipfile.ZipFile(zipname, 'w') as zipMe:        
            for f, file in enumerate(lista_files):
                if console_summary:
                    print('Zipping '+str(f+1)+' of '+str(nscraped)+': '+file)
                try:
                    zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)
                except (FileNotFoundError, OSError):
                    if console_summary:
                        print('\nSkipping '+file+', Not found\n' )
            
        if console_summary:
            print('\n'+'Finished - Zip File Created: '+zipname)
        if higher_dir:
            for file in lista_files:
                if '..' in file:
                    higher_files.append(file)
            if console_summary:
                print('\n********** WARNING: Files in parent directory(s) not zipped **********')
                for file in higher_files:
                    print(file)
                print('\nPlease add these manually to the zip file')

    if not console_summary:
        return missing
    else:
        return

def ensure_numpy_array(data):
    if isinstance(data, list):
        return np.array(data)
    elif isinstance(data, np.ndarray):
        return data
    else:
        raise ValueError("Input is neither a list nor a numpy array.")
                               
def rr_solver(
    zi: Union[npt.ArrayLike, list], ki: Union[npt.ArrayLike, list]
) -> Tuple[int, np.ndarray, np.ndarray, float, float]:
    """
    Solves for the root of the Rachford-Rice equation using a method that 
    gracefully handles catastrophic roundoff errors.
    The method is outlined in 2022 'Fluid Phase Equilibria' paper by M. Nielsen & H. Lia

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
    zi = ensure_numpy_array(zi)
    ki = ensure_numpy_array(ki)

    if len(zi) == 0 or len(ki) == 0:
        raise ValueError("zi and ki arrays must not be empty")
    if len(zi) != len(ki):
        raise ValueError(f"zi and ki arrays must have the same length (got {len(zi)} and {len(ki)})")

    zi_sum = np.sum(zi)
    if zi_sum <= 0:
        raise ValueError("Sum of zi must be positive")
    zi = zi / zi_sum  # Normalize feed compositions
    
    def rr(V: float) -> float:
        return np.dot(zi, (ki - 1) / (1 + V * (ki - 1)))
    
    def h(b: float) -> float:    # Eq 12b
        return np.sum(zi*b/(1+b*(phi_min - ci)))
    
    def dh(b: float) -> float:   # Eq 16b
        return np.sum(zi/(1 + b*(phi_min - ci)) ** 2)
    
    N_it = 0
    
    # Assess if solution nearest vapor or liquid solutions
    # and appropriately assign pseudo phase and ki's
    near_vapor = False
    if rr(0.5) > 0:
        near_vapor = True
    
    ki_hat = ki
    if near_vapor:
        ki_hat = 1/ki
    
    ci = np.where(np.abs(1 - ki_hat) > 1e-10, 1/(1-ki_hat), 0.0)  # Eq 10
        
    # Will use the transformed solution in terms of b (Eq 14b)
    phi_max = min(1/(1-min(ki_hat)), 0.5)    # Eq 11a
    phi_min = 1/(1-max(ki_hat))              # Eq 11b with upper 0.5 limit given we have chosen this phase
    b_min = 1/(phi_max - phi_min)            # Eq 15
    b_max = np.inf
    
    b = 1/(0.25 - phi_min)
    N_it = 0
    h_b = np.inf

    while abs(h_b) > EPS_T:
        N_it += 1

        h_b = h(b)
        dh_b = dh(b)

        if h_b > 0:
            b_max = b
        else:
            b_min = b

        b = b - h_b / dh_b

        if b < b_min or b > b_max:
            b = (b_min + b_max) / 2

        if N_it > MAX_ITR:
            break
    
    ui = -zi*ci*b/(1+b*(phi_min-ci))    # Eq 27b 
    phi = (1+b*phi_min)/b               # Rearranged Eq 14b
    
    if near_vapor:
        L = phi
        V = 1 - L
        yi = ui
        xi = ki_hat * ui                # Eq 28
        
    else:
        V = phi
        L = 1 - V
        xi = ui
        yi = ki_hat * ui                # Eq 28
        
    return N_it, yi, xi, V, L


# =============================================================================
# VFP Table Generation
# =============================================================================

def _format_vfpinj(table_num, datum_depth, flo_type, flo_rates, thp_values, bhp_array):
    """Format Eclipse VFPINJ keyword string."""
    lines = []
    lines.append("-- Generated by pyResToolbox make_vfpinj")
    lines.append("VFPINJ")
    lines.append(f"  {table_num}  {datum_depth:.5E}  {flo_type}  THP  FIELD  BHP /")
    lines.append("  " + "  ".join(f"{r:.5E}" for r in flo_rates) + " /")
    lines.append("  " + "  ".join(f"{t:.5E}" for t in thp_values) + " /")
    for it in range(len(thp_values)):
        bhps = bhp_array[it, :]
        bhp_str = "  ".join(f"{b:.5E}" for b in bhps)
        lines.append(f"  {it+1}  {bhp_str} /")
    return "\n".join(lines) + "\n"


def _format_vfpprod(table_num, datum_depth, flo_type, wfr_type, gfr_type,
                    flo_rates, thp_values, wfr_values, gfr_values, alq_values,
                    bhp_array):
    """Format Eclipse VFPPROD keyword string."""
    lines = []
    lines.append("-- Generated by pyResToolbox make_vfpprod")
    lines.append("VFPPROD")
    lines.append(f"  {table_num}  {datum_depth:.5E}  {flo_type}  {wfr_type}  {gfr_type}"
                 f"  THP  ''  FIELD  BHP /")
    lines.append("  " + "  ".join(f"{r:.5E}" for r in flo_rates) + " /")
    lines.append("  " + "  ".join(f"{t:.5E}" for t in thp_values) + " /")
    lines.append("  " + "  ".join(f"{w:.5E}" for w in wfr_values) + " /")
    lines.append("  " + "  ".join(f"{g:.5E}" for g in gfr_values) + " /")
    lines.append("  " + "  ".join(f"{a:.5E}" for a in alq_values) + " /")
    # BHP records: iterate ALQ(outer) > GFR > WFR > THP(inner)
    nthp = len(thp_values)
    nwfr = len(wfr_values)
    ngfr = len(gfr_values)
    nalq = len(alq_values)
    for ia in range(nalq):
        for ig in range(ngfr):
            for iw in range(nwfr):
                for it in range(nthp):
                    bhps = bhp_array[it, iw, ig, ia, :]
                    bhp_str = "  ".join(f"{b:.5E}" for b in bhps)
                    lines.append(f"  {it+1}  {iw+1}  {ig+1}  {ia+1}  {bhp_str} /")
    return "\n".join(lines) + "\n"


def make_vfpinj(
    table_num: int,
    completion,
    flo_type: str = 'WAT',
    vlpmethod: str = 'HB',
    flo_rates: list = None,
    thp_values: list = None,
    gas_pvt=None,
    gsg: float = 0.65,
    wsg: float = 1.07,
    oil_pvt=None,
    pb: float = 0,
    rsb: float = 0,
    sgsp: float = 0.65,
    api: float = 35,
    datum_depth: float = 0,
    export: bool = False,
    filename: str = '',
) -> dict:
    """ Generates an Eclipse VFPINJ keyword table for injection wells.

        Computes BHP as a function of flow rate and tubing head pressure using
        the nodal module VLP correlations.

        table_num: Integer VFP table number (>= 1)
        completion: A nodal.Completion object describing the wellbore
        flo_type: Flow rate type - 'WAT' (water), 'GAS' (gas), or 'OIL' (oil)
        vlpmethod: VLP correlation - 'HB', 'WG', 'GRAY', or 'BB'
        flo_rates: List of flow rates in ascending order.
                   Units: stb/d (WAT/OIL) or Mscf/d (GAS). Default provided if None
        thp_values: List of tubing head pressures (psia) in ascending order. Default provided if None
        gas_pvt: Optional GasPVT object for gas injection
        gsg: Gas specific gravity. Default 0.65
        wsg: Water specific gravity. Default 1.07
        oil_pvt: Optional OilPVT object for oil injection
        pb: Bubble point pressure (psia). For oil injection
        rsb: Solution GOR at Pb (scf/stb). For oil injection
        sgsp: Separator gas SG. Default 0.65
        api: API gravity. Default 35
        datum_depth: Bottom hole datum depth (ft). Default = completion total TVD
        export: If True, writes an Eclipse include file. Default False
        filename: Custom output filename. Default: VFPINJ_{table_num}.INC

        Returns: Dictionary with keys:
            table_num, datum_depth, flo_type, flo_rates, thp_values,
            bhp (2D numpy array shape NTHP x NFLO),
            n_failed, eclipse_string
    """
    from pyrestoolbox.nodal import fbhp as nodal_fbhp

    # Validation
    if not isinstance(table_num, int) or table_num < 1:
        raise ValueError(f"table_num must be a positive integer, got {table_num}")

    flo_type = flo_type.upper()
    if flo_type not in ('WAT', 'GAS', 'OIL'):
        raise ValueError(f"flo_type must be 'WAT', 'GAS', or 'OIL', got '{flo_type}'")

    vlpm = validate_methods(["vlpmethod"], [vlpmethod])

    if datum_depth <= 0:
        datum_depth = completion.total_tvd

    # Defaults
    if flo_rates is None:
        if flo_type == 'GAS':
            flo_rates = [1000, 5000, 10000, 20000, 50000]
        else:
            flo_rates = [100, 500, 1000, 2000, 5000, 10000]
    if thp_values is None:
        thp_values = [100, 200, 500, 1000, 2000]

    flo_rates = list(flo_rates)
    thp_values = list(thp_values)

    if len(flo_rates) < 2:
        raise ValueError("flo_rates must contain at least 2 values")
    if len(thp_values) < 1:
        raise ValueError("thp_values must contain at least 1 value")

    # Compute BHP grid
    nflo = len(flo_rates)
    nthp = len(thp_values)
    bhp_array = np.zeros((nthp, nflo))
    n_failed = 0

    for it, thp in enumerate(thp_values):
        for iflo, flo in enumerate(flo_rates):
            try:
                if flo_type == 'GAS':
                    bhp_val = nodal_fbhp(
                        thp=thp, completion=completion,
                        vlpmethod=vlpm, well_type='gas',
                        gas_pvt=gas_pvt,
                        qg_mmscfd=flo / 1000.0, cgr=0, qw_bwpd=0,
                        wsg=wsg, gsg=gsg, injection=True)
                elif flo_type == 'WAT':
                    bhp_val = nodal_fbhp(
                        thp=thp, completion=completion,
                        vlpmethod=vlpm, well_type='oil',
                        qt_stbpd=flo, gor=0, wc=1.0,
                        wsg=wsg, gsg=gsg, injection=True,
                        pb=0, rsb=0, sgsp=sgsp, api=api)
                else:  # OIL
                    bhp_val = nodal_fbhp(
                        thp=thp, completion=completion,
                        vlpmethod=vlpm, well_type='oil',
                        oil_pvt=oil_pvt,
                        qt_stbpd=flo, gor=0, wc=0,
                        wsg=wsg, gsg=gsg, injection=True,
                        pb=pb, rsb=rsb, sgsp=sgsp, api=api)

                bhp_array[it, iflo] = bhp_val
            except (RuntimeError, ValueError, ZeroDivisionError):
                bhp_array[it, iflo] = 1e-6
                n_failed += 1

    if n_failed > 0:
        print(f"Warning: {n_failed} of {nthp * nflo} VFPINJ BHP calculations failed")

    eclipse_str = _format_vfpinj(table_num, datum_depth, flo_type,
                                  flo_rates, thp_values, bhp_array)

    if export:
        fname = filename if filename else f"VFPINJ_{table_num}.INC"
        with open(fname, 'w') as f:
            f.write(eclipse_str)

    return {
        "table_num": table_num,
        "datum_depth": datum_depth,
        "flo_type": flo_type,
        "flo_rates": flo_rates,
        "thp_values": thp_values,
        "bhp": bhp_array,
        "n_failed": n_failed,
        "eclipse_string": eclipse_str,
    }


def make_vfpprod(
    table_num: int,
    completion,
    well_type: str = 'gas',
    vlpmethod: str = 'HB',
    flo_rates: list = None,
    thp_values: list = None,
    wfr_values: list = None,
    gfr_values: list = None,
    alq_values: list = None,
    gas_pvt=None,
    gsg: float = 0.65,
    oil_vis: float = 1.0,
    api: float = 45,
    pr: float = 0,
    oil_pvt=None,
    pb: float = 0,
    rsb: float = 0,
    sgsp: float = 0.65,
    wsg: float = 1.07,
    datum_depth: float = 0,
    export: bool = False,
    filename: str = '',
) -> dict:
    """ Generates an Eclipse VFPPROD keyword table for production wells.

        Computes BHP as a function of flow rate, tubing head pressure,
        water fraction, gas fraction, and artificial lift quantity using
        the nodal module VLP correlations.

        For gas wells: FLO=GAS (Mscf/d), WFR=WGR (stb/Mscf), GFR=OGR (stb/Mscf)
        For oil wells: FLO=OIL (stb/d), WFR=WCT (fraction), GFR=GOR (Mscf/stb)

        table_num: Integer VFP table number (>= 1)
        completion: A nodal.Completion object describing the wellbore
        well_type: 'gas' or 'oil'
        vlpmethod: VLP correlation - 'HB', 'WG', 'GRAY', or 'BB'
        flo_rates: List of flow rates in ascending order. Default provided if None
        thp_values: List of tubing head pressures (psia). Default provided if None
        wfr_values: List of water fraction values. Default provided if None
                    Gas wells: WGR in stb/Mscf. Oil wells: WCT as fraction (0-1)
        gfr_values: List of gas fraction values. Default provided if None
                    Gas wells: OGR in stb/Mscf. Oil wells: GOR in Mscf/stb
        alq_values: List of artificial lift values. Default [0]
        gas_pvt: Optional GasPVT object
        gsg: Gas specific gravity. Default 0.65
        oil_vis: Condensate viscosity (cP). Default 1.0
        api: API gravity (condensate for gas wells, oil for oil wells). Default 45
        pr: Reservoir pressure (psia) for condensate dropout model. 0 disables. Default 0
        oil_pvt: Optional OilPVT object for oil wells
        pb: Bubble point pressure (psia). Required for oil wells
        rsb: Solution GOR at Pb (scf/stb). Required for oil wells
        sgsp: Separator gas SG. Default 0.65
        wsg: Water specific gravity. Default 1.07
        datum_depth: Bottom hole datum depth (ft). Default = completion total TVD
        export: If True, writes an Eclipse include file. Default False
        filename: Custom output filename. Default: VFPPROD_{table_num}.INC

        Returns: Dictionary with keys:
            table_num, datum_depth, well_type, flo_type, wfr_type, gfr_type,
            flo_rates, thp_values, wfr_values, gfr_values, alq_values,
            bhp (5D numpy array shape NTHP x NWFR x NGFR x NALQ x NFLO),
            n_failed, eclipse_string
    """
    from pyrestoolbox.nodal import fbhp as nodal_fbhp

    # Validation
    if not isinstance(table_num, int) or table_num < 1:
        raise ValueError(f"table_num must be a positive integer, got {table_num}")

    well_type = well_type.lower()
    if well_type not in ('gas', 'oil'):
        raise ValueError(f"well_type must be 'gas' or 'oil', got '{well_type}'")

    vlpm = validate_methods(["vlpmethod"], [vlpmethod])

    if datum_depth <= 0:
        datum_depth = completion.total_tvd

    # Defaults
    if flo_rates is None:
        flo_rates = [1000, 5000, 10000, 20000, 50000, 100000] if well_type == 'gas' \
                    else [100, 500, 1000, 2000, 5000, 10000]
    if thp_values is None:
        thp_values = [100, 200, 500, 1000, 2000]
    if wfr_values is None:
        wfr_values = [0, 1, 5, 10] if well_type == 'gas' else [0, 0.2, 0.5, 0.8]
    if gfr_values is None:
        gfr_values = [0, 10, 50, 100] if well_type == 'gas' else [0.2, 0.5, 1.0, 2.0]
    if alq_values is None:
        alq_values = [0]

    flo_rates = list(flo_rates)
    thp_values = list(thp_values)
    wfr_values = list(wfr_values)
    gfr_values = list(gfr_values)
    alq_values = list(alq_values)

    if len(flo_rates) < 2:
        raise ValueError("flo_rates must contain at least 2 values")

    # Validate WCT range for oil wells
    if well_type == 'oil':
        if any(w >= 1.0 for w in wfr_values):
            raise ValueError("WCT values must be less than 1.0 for oil production wells")

    # Type labels
    if well_type == 'gas':
        flo_type, wfr_type, gfr_type = 'GAS', 'WGR', 'OGR'
    else:
        flo_type, wfr_type, gfr_type = 'OIL', 'WCT', 'GOR'

    # Compute BHP grid
    nflo = len(flo_rates)
    nthp = len(thp_values)
    nwfr = len(wfr_values)
    ngfr = len(gfr_values)
    nalq = len(alq_values)
    total = nflo * nthp * nwfr * ngfr * nalq

    bhp_array = np.zeros((nthp, nwfr, ngfr, nalq, nflo))
    n_failed = 0
    n_computed = 0

    for ia, alq in enumerate(alq_values):
        for ig, gfr in enumerate(gfr_values):
            for iw, wfr in enumerate(wfr_values):
                for it, thp in enumerate(thp_values):
                    for iflo, flo in enumerate(flo_rates):
                        n_computed += 1
                        try:
                            if well_type == 'gas':
                                qg_mmscfd = flo / 1000.0
                                qw_bwpd = wfr * flo  # WGR(stb/Mscf) * rate(Mscf/d)
                                cgr = gfr * 1000.0    # OGR(stb/Mscf) -> stb/MMscf
                                bhp_val = nodal_fbhp(
                                    thp=thp, completion=completion,
                                    vlpmethod=vlpm, well_type='gas',
                                    gas_pvt=gas_pvt,
                                    qg_mmscfd=qg_mmscfd, cgr=cgr,
                                    qw_bwpd=qw_bwpd, oil_vis=oil_vis,
                                    api=api, pr=pr, wsg=wsg, gsg=gsg)
                            else:
                                wct = wfr
                                gor_scf_stb = gfr * 1000.0
                                qt_stbpd = flo / (1.0 - wct) if wct > 0 else flo
                                bhp_val = nodal_fbhp(
                                    thp=thp, completion=completion,
                                    vlpmethod=vlpm, well_type='oil',
                                    oil_pvt=oil_pvt,
                                    qt_stbpd=qt_stbpd, gor=gor_scf_stb,
                                    wc=wct, wsg=wsg, gsg=gsg,
                                    pb=pb, rsb=rsb, sgsp=sgsp, api=api)

                            bhp_array[it, iw, ig, ia, iflo] = bhp_val
                        except (RuntimeError, ValueError, ZeroDivisionError):
                            bhp_array[it, iw, ig, ia, iflo] = 1e-6
                            n_failed += 1

    if n_failed > 0:
        print(f"Warning: {n_failed} of {total} VFPPROD BHP calculations failed")

    eclipse_str = _format_vfpprod(table_num, datum_depth, flo_type, wfr_type, gfr_type,
                                   flo_rates, thp_values, wfr_values, gfr_values,
                                   alq_values, bhp_array)

    if export:
        fname = filename if filename else f"VFPPROD_{table_num}.INC"
        with open(fname, 'w') as f:
            f.write(eclipse_str)

    return {
        "table_num": table_num,
        "datum_depth": datum_depth,
        "well_type": well_type,
        "flo_type": flo_type,
        "wfr_type": wfr_type,
        "gfr_type": gfr_type,
        "flo_rates": flo_rates,
        "thp_values": thp_values,
        "wfr_values": wfr_values,
        "gfr_values": gfr_values,
        "alq_values": alq_values,
        "bhp": bhp_array,
        "n_failed": n_failed,
        "eclipse_string": eclipse_str,
    }


# ============================================================================
# Black Oil Table Generation (entry point in simtools, logic in oil module)
# ============================================================================

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
    pbmethod: pb_method = pb_method.VELAR,
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
    wt: Salt wt% (0-100) in brine. Default = 0
    ch4_sat: Degree of methane saturation (0 - 1) in brine. Default = 0
    export: Boolean flag that controls whether to export full table to excel, and separate PVDG and PVDO include files. Default is False
    pvto: Boolean flag that controls whether the pvto live oil Eclipse format will be generated. This can only be active if export flag is also True;
          - extends bubble point line up to maximum pressure
          - generates undersaturated oil propeties
          - writes out PVTO include file
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

    sg_g, sg_sp = oil.check_sgs(sg_g=sg_g, sg_sp=0)
    pmin = max(pmin, psc)
    sg_o = oil.oil_sg(api)

    # Step 1: Resolve Pb, Rsb, and scaling factor
    pb, rsb, rsb_frac, rsb_max, pb_i, rsb_i = _oil_impl._resolve_pb_rsb(
        pb, rsb, degf, api, sg_sp, sg_g, pvto, pmax, rsmethod, pbmethod
    )

    # Step 2: Build pressure array
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

    # Step 3: Compute all PVT properties
    (rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
     usat_p, usat_bo, usat_uo) = _oil_impl._build_bot_tables(
        pressures, pb, rsb, rsb_frac, rsb_max, sg_o, sg_g, sg_sp,
        api, degf, pvto, wt, ch4_sat,
        zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod,
    )

    # Step 4: Assemble results and optionally export
    return _oil_impl._format_bot_results(
        pressures, rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
        usat_p, usat_bo, usat_uo, sg_o, sg_g, pi, degf, wt, ch4_sat,
        pb_i, rsb_i, rsb_frac, pvto, export, zmethod, cmethod,
    )


# ============================================================================
# PVTW Table Generation (moved from brine module for v3.0)
# ============================================================================

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
        bw, lden, visw, cw, rsw = brine.brine_props(p=p, degf=degf, wt=wt, ch4_sat=ch4_sat)
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
    bw_ref, den_ref, visw_ref, cw_ref, rsw_ref = brine.brine_props(
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