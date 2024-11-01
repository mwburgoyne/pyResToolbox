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
import os
from os.path import exists
import zipfile

from typing import Union, List, Tuple
import numpy as np
import numpy.typing as npt

import pandas as pd
from tabulate import tabulate
from gwr_inversion import gwr
from mpmath import mp

EPS_T = 1e-15
MAX_ITR = 100

class kr_family(Enum):  # Relative permeability family type
    COR = 0
    LET = 1


class kr_table(Enum):  # Relative permeability table type
    SWOF = 0
    SGOF = 1
    SGWFN = 2
    

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
            print("File name needs to be an IX print file with .PRT extension")
            return

    if filename == "":  # Show selection in local directory
        prt_files = glob.glob("*.PRT", recursive=False)
        if len(prt_files) == 0:
            print("No .PRT files exist in this directory - Terminating script")
            sys.exit()

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
                print(
                    "\nIndex entered outside range permitted - Terminating script"
                )
                sys.exit()
        else:
            prt_file_idx = 0

        filename = prt_files[prt_file_idx]

    if not silent:
        print("Processing " + filename + "\n")
    file1 = open(filename, "r")
    count = 0
    grab_line1 = False
    grab_line2 = False
    max_it = 12
    timesteps = []
    tables = []

    while True:
        count += 1
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
        print("Does not appear to be a valid IX PRT file")
        return

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
        return occurence_count.most_common(1)[0][0]

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
            problem_data.append(most_frequent(dat))
            problem_data_count.append(dat.count(problem_data[-1]))
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
        except:
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
    return s ** L / (s ** L + E * (1 - s) ** T)
    
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
    export: bool = False,
) -> pd.DataFrame:
    """ Returns ECLIPSE styled relative permeability tables
        Users need only define parameters relevant to their table / family selection
        rows: Integer value specifying the number of table rows desired
        krtable: A string or kr_table Enum class that specifies one of three table type choices;
                   SWOF: Water / Oil table
                   SGOF: Gas / Oil table
                   SGFN: Gas / Water table
        krfamily: A string or kr_family Enum class that specifies one of two curve function choices;
                   COR: Corey Curve function
                   LET: LET Relative permeability function
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
        export: Boolean value that controls whether an include file with same name as krtable is created. Default: False
    """

    if type(krtable) == str:
        try:
            krtable = kr_table[krtable.upper()]
        except:
            print("Incorrect table type specified")
            sys.exit()
    if type(krfamily) == str:
        try:
            krfamily = kr_family[krfamily.upper()]
        except:
            print("Incorrect krfamily specified")
            sys.exit()

    def kr_SWOF(
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
    ) -> pd.DataFrame:

        if no * nw <= 0:  # Not enough information for Corey curves
            corey_info = False
        else:
            corey_info = True
        if (
            Lw * Ew * Tw * Lo * Eo * To <= 0
        ):  # Not enough information for LET curves
            let_info = False
        else:
            let_info = True

        ndiv = rows
        if swcr > swc:
            ndiv -= 2
        if sorw > 0:
            ndiv -= 1
        ndiv = min(ndiv, rows - 1)

        sw_eps = [swc, swcr, 1 - sorw, 1]
        swn = np.arange(0, 1, 1 / ndiv)
        sw = swn * (1 - swcr - sorw) + swcr
        sw = list(sw) + sw_eps
        sw = list(set(sw))
        sw.sort()
        sw = np.array(sw)

        # Assign water relative permeabilities
        swn = (sw - swcr) / (1 - swcr - sorw)
        swn = np.clip(swn, 0, 1)

        if krfamily.name == "COR":
            if not corey_info:
                print(
                    "Not enough information for SWOF Corey Curves. Check if no and nw are defined"
                )
                return
            krw = krwmax*corey(swn, nw)
            if sorw > 0:
                krw[-1] = 1
        if krfamily.name == "LET":
            if not let_info:
                print(
                    "Not enough information for SWOF LET Curves. Check if Lw, Ew, Tw, Lo, Eo & To are defined"
                )
                return
            krw = krwmax*LET(swn, Lw, Ew, Tw)
            if sorw > 0:
                krw[-1] = 1

        # Assign oil relative permeabilities
        swn = (sw - swc) / (1 - swc - sorw)
        swn = np.clip(swn, 0, 1)
        if krfamily.name == "COR":
            kro = kromax*corey(1 - swn, no)
        if krfamily.name == "LET":
            kro = kromax*LET(1 - swn, Lo, Eo, To)

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

    def kr_SGOF(
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
    ) -> pd.DataFrame:

        if not no * ng:  # Not enough information for Corey curves
            corey_info = False
        else:
            corey_info = True
        if (
            not Lg * Eg * Tg * Lo * Eo * To
        ):  # Not enough information for LET curves
            let_info = False
        else:
            let_info = True

        ndiv = rows
        if sgcr > 0:
            ndiv -= 2
        if sorg > 0:
            ndiv -= 1
        ndiv = min(ndiv, rows - 1)

        sg_eps = [0, 1 - swc - sorg]
        sgn = np.arange(0, 1 + 1 / ndiv, 1 / ndiv)
        sg = sgn * (1 - swc - sorg)
        sg = list(sg) + sg_eps
        sg = list(set(sg))
        sg.sort()
        sg = np.array(sg)

        # Assign gas relative permeabilities
        sgn = sg / (1 - swc - sorg)
        sgn = np.clip(sgn, 0, 1)
        if krfamily.name == "COR":
            if not corey_info:
                print(
                    "Not enough information for SGOF Corey Curves. Check if no and ng are defined"
                )
                return
            krg = krgmax * corey(sgn,ng)
        if krfamily.name == "LET":
            if not let_info:
                print(
                    "Not enough information for SGOF LET Curves. Check if Lg, Eg, Tg, Lo, Eo & To are defined"
                )
                return
            krg = krgmax*LET(sgn, Lg, Eg, Tg)

        # Assign oil relative permeabilities
        if krfamily.name == "COR":
            kro = kromax * corey(1 - sgn, no)
        if krfamily.name == "LET":
            kro = kromax*LET(1 - sgn, Lo, Eo, To)

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

    def kr_SGWFN(
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
    ) -> pd.DataFrame:
        if ng * nw <= 0:  # Not enough information for Corey curves
            corey_info = False
        else:
            corey_info = True
        if (
            Lw * Ew * Tw * Lg * Eg * Tg <= 0
        ):  # Not enough information for LET curves
            let_info = False
        else:
            let_info = True

        ndiv = rows
        if sgcr > 0:
            ndiv -= 1
        ndiv = min(ndiv, rows - 1)

        sg_eps = [0, sgcr, 1 - swc]
        ndiv = rows - 1
        sgn = np.arange(0, 1, 1 / ndiv)
        sg = sgn * (1 - swc - sgcr) + sgcr
        sg = list(sg) + sg_eps
        sg = list(set(sg))
        sg.sort()
        sg = np.array(sg)

        # Assign gas relative permeabilities
        sgn = (sg - sgcr) / (1 - swc - sgcr)
        sgn = np.clip(sgn, 0, 1)
        if krfamily.name == "COR":
            if not corey_info:
                print(
                    "Not enough information for SGWFN Corey Curves. Check if nw and ng are defined"
                )
                return
            krg = krgmax * corey(sgn, ng)
            
        if krfamily.name == "LET":
            if not let_info:
                print(
                    "Not enough information for SGWFN LET Curves. Check if Lg, Eg, Tg, Lw, Ew & Tw are defined"
                )
                return
            krg = krgmax*LET(sgn, Lg, Eg, Tg)

        # Assign water relative permeabilities
        sgn = (sg) / (1 - swc)
        sgn = np.clip(sgn, 0, 1)
        if krfamily.name == "COR":
            krw = krwmax * corey(1 - sgn, nw)
        if krfamily.name == "LET":
            krw = krwmax*LET(1-sgn, Lw, Ew, Tw)


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

    # Consistency checks
    fail = False
    swcr = max(swc, swcr)
    if sorg + sgcr + swc >= 1:
        print("sorg+sgcr+swc must be less than 1")
        fail = True
    if sorg + sgcr + swc >= 1:
        print("sorg+sgcr+swc must be less than 1")
        fail = True
    if sorw + swcr >= 1:
        print("sorw+swcr must be less than 1")
        fail = True
    if fail:
        print("Saturation consistency check failure: Check your inputs")
        sys.exit()

    if krtable.name == "SWOF":
        return kr_SWOF(
            rows,
            krtable,
            krfamily,
            kromax,
            krgmax,
            krwmax,
            swc,
            swcr,
            sorg,
            sorw,
            sgcr,
            no,
            nw,
            ng,
            Lw,
            Ew,
            Tw,
            Lo,
            Eo,
            To,
            Lg,
            Eg,
            Tg,
        )
    if krtable.name == "SGOF":
        return kr_SGOF(
            rows,
            krtable,
            krfamily,
            kromax,
            krgmax,
            krwmax,
            swc,
            swcr,
            sorg,
            sorw,
            sgcr,
            no,
            nw,
            ng,
            Lw,
            Ew,
            Tw,
            Lo,
            Eo,
            To,
            Lg,
            Eg,
            Tg,
        )
    if krtable.name == "SGWFN":
        return kr_SGWFN(
            rows,
            krtable,
            krfamily,
            kromax,
            krgmax,
            krwmax,
            swc,
            swcr,
            sorg,
            sorw,
            sgcr,
            no,
            nw,
            ng,
            Lw,
            Ew,
            Tw,
            Lo,
            Eo,
            To,
            Lg,
            Eg,
            Tg,
        )
    print("Check that you have specified table type as SWOF, SGOF or SGWFN")
    sys.exit()

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
    if min(ReDs) <=1:
        print('\n\n\n -------  ReDs must all be strictly greater than 1.0  -------\n\n\n')
        return 

    # Eq 16 from SPE 81428
    #def laplace_Qs(s: float, ReD: float):
    #    x = mp.sqrt(s)
    #    # pre-calculate duplicated bessel functions for greater efficiency
    #    i1sReD = mp.besseli(1, x * ReD)
    #    k1sReD = mp.besselk(1, x * ReD)
    #    numerator = i1sReD * mp.besselk(1, x) - k1sReD * mp.besseli(1, x)
    #    denominator = s ** 1.5 * (
    #        k1sReD * mp.besseli(0, x) + i1sReD * mp.besselk(0, x)
    #    )
    #    return numerator / denominator

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

        text_file = open("INFLUENCE.INC", "w")
        text_file.write(inc_out)
        text_file.close()

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
            print( 'No '+mask+' files exist in this directory - Terminating script')
            sys.exit()
    
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
            print(' ')
            print( 'Index entered outside range permitted - Terminating script')
            sys.exit()
        
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
    missing_parents = []
    # Start stepping through, looking for INCLUDE file statements
    get_include = False
    got_all = False
    
    nscraped = 0
    higher_dir = False
    missing = []
    endfound = False
    while not got_all:
        if console_summary:
            print('Scanning through: '+files2scrape[nscraped])
                    
        # Load file into list
        try:
            lines = list(Tuple(open(files2scrape[nscraped], 'r')))
        except:
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
                endfound = True
                continue
            if line.upper()[:7]=='INCLUDE':
                get_include = True

            if get_include:
                
                # Remove any trailing comments that might give false negatives
                line = line.split('--')[0]
                line = line.split('#')[0]
                
                line = line.replace('"',"'") # Remove double inverted commas, replace with single commas
                if '.' not in line and "'" not in line: # Filename not in this line
                    continue
                
                include_file = line.split("'")[1] # Get filename
                include_file=include_file.strip()
                if include_file not in files2scrape:
                    files2scrape.append(include_file)   
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
                    sys.exit()
        lista_files = files2scrape
        dotindex = ''.join(lista_files[0]).rindex('.')
        zipname = lista_files[0][:dotindex]+'.zip'
        with zipfile.ZipFile(zipname, 'w') as zipMe:        
            for f, file in enumerate(lista_files):
                if console_summary:
                    print('Zipping '+str(f+1)+' of '+str(nscraped)+': '+file)
                try:
                    zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)
                except:
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
    
    zi = zi/np.sum(zi) # Normalize feed compositions
    
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
    
    ci = 1/(1-ki_hat)  # Eq 10
        
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