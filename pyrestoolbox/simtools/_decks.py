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

Simulation deck utilities: IX PRT convergence parsing and recursive
ECL/IX INCLUDE file checking/zipping.
"""

from collections import Counter
from dataclasses import dataclass
import glob
import zipfile
from os.path import exists

import pandas as pd
from tabulate import tabulate


def ix_extract_problem_cells(filename: str = "", silent: bool = False, non_interactive: bool = True) -> list:
    """
    Processes Intersect PRT file to extract convergence issue information
    Prints a summary of worst offenders to terminal (if silent=False), and returns a list
    of sorted dataframes summarising all entities in final convergence row in the PRT files
    List returned is [well_pressure_df, grid_pressure_df, sat_change_df, comp_change_df]
    filename: If empty, will search local directory for PRT file and present list to select from if more than one exists.
              If a filename is furnished, or only one file exists, then no selection will be presented
    silent: False will return only the list of dataframes, with nothing echoed to the terminal
            True will return summary of worst entities to the terminal
    non_interactive: If True (default), raise ValueError instead of prompting via input() when multiple .PRT files are found
                     and no filename was supplied. Suits scripts/agents where stdin is not available.
                     Pass non_interactive=False for the legacy interactive prompt.
    """

    if filename != "":  # A Filename has been provided
        if "PRT" not in filename.upper():
            raise ValueError("File name needs to be an IX print file with .PRT extension: " + filename)

    if filename == "":  # Show selection in local directory
        prt_files = glob.glob("*.PRT", recursive=False)
        if len(prt_files) == 0:
            raise FileNotFoundError("No .PRT files exist in this directory")

        if len(prt_files) > 1:
            if non_interactive:
                raise ValueError(
                    f"Multiple .PRT files found in directory ({prt_files}); "
                    "pass filename=... explicitly in non-interactive mode"
                )
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


@dataclass
class _DeckScan:
    """Result of crawling a deck's INCLUDE tree."""
    files: list          # every file reached, deck first
    parents: list        # the file each entry was referenced from
    missing: list        # entries that do not exist on disk
    missing_parents: list
    higher_dir: bool     # any INCLUDE resolved above the deck directory


def _select_deck_files(mask):
    """Prompt the user to pick deck files matching mask from the current directory."""
    input_files = []
    for files in mask:
        for file in glob.glob(files.upper()):  # Grab all files in the directory with this mask
            input_files.append(file)
        for file in glob.glob(files.lower()):  # In case lowercase extension used
            input_files.append(file)
    input_files = list(set(input_files))  # Upper/lower check can yield duplicates

    if len(input_files) == 0:
        raise FileNotFoundError('No '+', '.join(mask)+' files exist in this directory')

    print(' ')
    input_files.sort()

    table = [[i, name] for i, name in enumerate(input_files)]
    print(tabulate(table, headers=['Index', 'File Name']))
    print(' ')
    file_idx = input('Please choose index(s) of file to parse separated by commas (0 - '+str(len(input_files)-1)+') :')
    file_idxs = [int(x) for x in file_idx.split(',')]

    if not all(0 <= item < len(input_files) for item in file_idxs):
        raise ValueError('Index entered outside range permitted')

    return [input_files[x] for x in file_idxs]


def _include_target(line):
    """Filename named on an INCLUDE continuation line, or None.

    Returns the quoted filename, or None when this line carries no name yet
    (a bare INCLUDE keyword, so the name is on a following line). Raises nothing:
    a malformed line is reported by the caller resetting its INCLUDE state.
    """
    line = line.split('--')[0]        # Drop trailing comments that give false negatives
    line = line.split('#')[0]
    line = line.replace('"', "'")     # Normalise double quotes to single
    if '.' not in line and "'" not in line:
        return None                   # Filename is not on this line
    parts = line.split("'")
    if len(parts) < 2:                # Malformed INCLUDE line
        return ''
    return parts[1].strip()


def _crawl_deck_includes(files2scrape, console_summary=True):
    """Walk a deck's INCLUDE tree breadth-first, to unlimited depth.

    Each newly named file is appended to the work list, so nested INCLUDEs are
    picked up without recursion. Files that cannot be opened are recorded as
    missing along with the file that referenced them.
    """
    files = list(files2scrape)
    parents = list(files2scrape)
    seen = set(files2scrape)
    missing = []
    missing_parents = []
    higher_dir = False

    nscraped = 0
    while nscraped < len(files):
        current = files[nscraped]
        if console_summary:
            print('Scanning through: '+current)

        try:
            lines = list(open(current, 'r'))
        except (FileNotFoundError, PermissionError, OSError):
            if not exists(current) and current not in missing:
                missing.append(current)
                missing_parents.append(parents[nscraped])
            nscraped += 1
            continue

        get_include = False
        for line in lines:
            line = line.strip()
            if line[:3] == '--' or line[:1] == '#':   # Skip all comments
                continue
            if line[:4] == 'END':                     # Skip everything after END
                break
            if line.upper()[:7] == 'INCLUDE':
                get_include = True
            if not get_include:
                continue

            include_file = _include_target(line)
            if include_file is None:                  # Name is on a later line
                continue
            get_include = False
            if not include_file or include_file in seen:
                continue
            files.append(include_file)
            seen.add(include_file)
            parents.append(current)
            if '..' in include_file:
                higher_dir = True

        nscraped += 1

    return _DeckScan(files=files, parents=parents, missing=missing,
                     missing_parents=missing_parents, higher_dir=higher_dir)


def _warn_higher_dir(files, heading, advice, console_summary):
    """Report INCLUDE files that sit above the deck directory."""
    higher = [f for f in files if '..' in f]
    if not higher or not console_summary:
        return
    print(heading)
    for file in higher:
        print(file)
    print(advice)


def _zip_deck_files(scan, console_summary, non_interactive):
    """Zip every scraped file alongside the deck, named after the deck."""
    if console_summary:
        print('\n'+str(len(scan.files))+' files to zip')
    if len(scan.missing) > 0:
        if non_interactive:
            raise RuntimeError(
                f"Cannot zip: {len(scan.missing)} INCLUDE file(s) missing: {scan.missing}"
            )
        if console_summary:
            cont = input('Continue to zip even with missing files? (Y/n): ')
            if cont.upper() == 'N':
                raise RuntimeError("Zip aborted by user due to missing files")

    dotindex = ''.join(scan.files[0]).rindex('.')
    zipname = scan.files[0][:dotindex]+'.zip'
    with zipfile.ZipFile(zipname, 'w') as zipMe:
        for f, file in enumerate(scan.files):
            if console_summary:
                print('Zipping '+str(f+1)+' of '+str(len(scan.files))+': '+file)
            try:
                zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)
            except (FileNotFoundError, OSError):
                if console_summary:
                    print('\nSkipping '+file+', Not found\n')

    if console_summary:
        print('\n'+'Finished - Zip File Created: '+zipname)
    if scan.higher_dir:
        _warn_higher_dir(
            scan.files,
            '\n********** WARNING: Files in parent directory(s) not zipped **********',
            '\nPlease add these manually to the zip file', console_summary)


def zip_check_sim_deck(files2scrape = [], tozip = True, console_summary = True, non_interactive = True):
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
        non_interactive: If True (default), raise ValueError instead of prompting via input(). Suits scripts/agents where stdin is not available.
                         In this mode: files2scrape must be provided, tozip is respected as given, and zipping with missing files raises rather than asking.
                         Pass non_interactive=False for the legacy interactive prompts.
    """

    types = ('*.DATA', '*.afi')
    if len(files2scrape) == 0:
        if non_interactive:
            raise ValueError(
                "files2scrape must be provided in non-interactive mode "
                "(no input() prompt available)"
            )
        files2scrape = _select_deck_files(types)
        tozip = True
        method = input('Zip or Check files? (Z/c): ')
        if method.upper() == 'C':
            tozip = False

    scan = _crawl_deck_includes(files2scrape, console_summary)

    if len(scan.missing) > 0 and console_summary:
        print('\n****** MISSING FILES ******\n')
        for f, file in enumerate(scan.missing):
            print(file, 'from', scan.missing_parents[f])

    if tozip:
        _zip_deck_files(scan, console_summary, non_interactive)
    elif len(scan.missing) == 0:
        if console_summary:
            print('\nALL INCLUDE FILES FOUND\n')
        _warn_higher_dir(
            scan.files,
            '\n********** WARNING: Some INCLUDE files in a parent directory **********',
            '\nThese would need to be added manually to a zip file', console_summary)

    if not console_summary:
        return scan.missing
    return
