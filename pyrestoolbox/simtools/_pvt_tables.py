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

PVT table generation: black oil tables (PVDO/PVDG/PVTO) and water PVT
tables (PVTW) in FIELD or Eclipse METRIC units.
"""

import numpy as np
import pandas as pd

from pyrestoolbox.classes import (z_method, c_method, pb_method, rs_method,
                                  bo_method, deno_method, co_method)
from pyrestoolbox.constants import (
    psc,
    CUFTperBBL,
    BAR_TO_PSI, PSI_TO_BAR,
    degc_to_degf,
    SCF_PER_STB_TO_SM3_PER_SM3, SM3_PER_SM3_TO_SCF_PER_STB,
    LBCUFT_TO_KGM3,
    INVPSI_TO_INVBAR,
)
from pyrestoolbox.validate import validate_methods
import pyrestoolbox.brine as brine
import pyrestoolbox.oil as oil
from pyrestoolbox.oil import _tables as _oil_tables

# Bg unit conversion: rb/mscf -> rm3/sm3 is a factor of 1000/5.614583 = 178.108
# (1000 scf/mscf divided by 5.614583 cuft/bbl)
_RB_PER_MSCF_TO_RM3_PER_SM3 = 1.0 / (1000.0 / CUFTperBBL)


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
    vis_frac: float = 1.0,
    metric: bool = False,
) -> dict:
    """
    Creates data required for Oil-Gas-Water black oil tables
    Returns dictionary of results, with index:
      - bot: Pandas table of blackoil data (for PVTO == False), or Saturated properties to pmax (if PVTO == True)
      - deno: ST Oil Density (lb/cuft or kg/m3 if metric)
      - deng: ST Gas Density (lb/cuft or kg/m3 if metric)
      - denw: Water Density at Pi (lb/cuft or kg/m3 if metric),
      - cw: Water Compressibility at Pi (1/psi or 1/bar if metric)
      - uw: Water Viscosity at Pi (cP))
      - pb: Bubble point pressure either calculated (if only Rsb provided), or supplied by user
      - rsb: Solution GOR at Pb either calculated (if only Pb provided), or supplied by user
      - rsb_scale: The scaling factor that was needed to match user supplied Pb and Rsb
      - usat: a list of understaurated values (if PVTO == True) [usat_p, usat_bo, usat_uo]. This will be empty if PVTO == False

    If user species Pb or Rsb only, the corresponding property will be calculated
    If both Pb and Rsb are specified, then Pb calculations will be adjusted to honor both

    pi: Initial reservoir pressure (psia / barsa if metric). Used to return water properties at initial pressure
    pb: Bubble point pressure (psia / barsa if metric)
    rsb: Oil solution GOR at Pb (scf/stb / sm3/sm3 if metric)
    degf: Reservoir Temperature (deg F / deg C if metric)
    sg_g: Weighted average specific gravity of surface gas (relative to air).
    api: Stock tank oil density (deg API).
    pmax: Maximum pressure to calcuate table to (psia / barsa if metric)
    pmin: Minimum pressure to calculate table to (psia / barsa if metric). Default = 25
    nrows: Number of BOT rows. Default = 20
    wt: Salt wt% (0-100) in brine. Default = 0
    ch4_sat: Degree of methane saturation (0 - 1) in brine. Default = 0
    export: Boolean flag that controls whether to export full table to excel, and separate PVDG and PVDO include files. Default is False
    pvto: Boolean flag that controls whether the pvto live oil Eclipse format will be generated. This can only be active if export flag is also True;
          - extends bubble point line up to maximum pressure
          - generates undersaturated oil propeties
          - writes out PVTO include file
    metric: If True, inputs/outputs use Eclipse METRIC units (barsa, degC, sm3/sm3, kg/m3).
            If False (default), inputs/outputs use FIELD units (psia, degF, scf/stb, lb/cuft).
    """
    # --- Metric input conversion ---
    if metric:
        pi = pi * BAR_TO_PSI
        pmax = pmax * BAR_TO_PSI
        pmin = pmin * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

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
    pb, rsb, rsb_frac, rsb_max, pb_i, rsb_i = _oil_tables._resolve_pb_rsb(
        pb, rsb, degf, api, sg_sp, sg_g, pvto, pmax, rsmethod, pbmethod
    )

    # Step 2: Build pressure array
    pmax = max(pb, pmax)
    pbi = pb
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
     usat_p, usat_bo, usat_uo) = _oil_tables._build_bot_tables(
        pressures, pb, rsb, rsb_frac, rsb_max, sg_o, sg_g, sg_sp,
        api, degf, pvto, wt, ch4_sat,
        zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod,
        vis_frac=vis_frac,
    )

    # Step 4: Assemble results and optionally export
    results = _oil_tables._format_bot_results(
        pressures, rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
        usat_p, usat_bo, usat_uo, sg_o, sg_g, pi, degf, wt, ch4_sat,
        pb_i, rsb_i, rsb_frac, pvto, export, zmethod, cmethod,
        vis_frac=vis_frac,
    )

    # --- Metric output conversion ---
    if metric:
        df = results["bot"]
        df["Pressure (psia)"] = df["Pressure (psia)"] * PSI_TO_BAR
        df.rename(columns={"Pressure (psia)": "Pressure (barsa)"}, inplace=True)
        df["Rs (mscf/stb)"] = df["Rs (mscf/stb)"] * 1000 * SCF_PER_STB_TO_SM3_PER_SM3  # mscf/stb -> scf/stb -> sm3/sm3
        df.rename(columns={"Rs (mscf/stb)": "Rs (sm3/sm3)"}, inplace=True)
        # Bo, Bg, Bw are dimensionless FVF ratios — no conversion needed
        df.rename(columns={"Bo (rb/stb)": "Bo (rm3/sm3)"}, inplace=True)
        df["Deno (lb/cuft)"] = df["Deno (lb/cuft)"] * LBCUFT_TO_KGM3
        df.rename(columns={"Deno (lb/cuft)": "Deno (kg/m3)"}, inplace=True)
        # uo, ug, uw in cP — no conversion
        df["Co (1/psi)"] = df["Co (1/psi)"] * INVPSI_TO_INVBAR
        df.rename(columns={"Co (1/psi)": "Co (1/bar)"}, inplace=True)
        df["Bg (rb/mscf"] = df["Bg (rb/mscf"] * _RB_PER_MSCF_TO_RM3_PER_SM3
        df.rename(columns={"Bg (rb/mscf": "Bg (rm3/sm3)"}, inplace=True)
        df["Cg (1/psi)"] = df["Cg (1/psi)"] * INVPSI_TO_INVBAR
        df.rename(columns={"Cg (1/psi)": "Cg (1/bar)"}, inplace=True)
        df.rename(columns={"Bw (rb/stb)": "Bw (rm3/sm3)"}, inplace=True)
        results["bot"] = df

        # Scalar outputs
        results["deno"] = results["deno"] * LBCUFT_TO_KGM3  # lb/cuft -> kg/m3
        results["deng"] = results["deng"] * LBCUFT_TO_KGM3
        results["denw"] = results["denw"] * LBCUFT_TO_KGM3
        results["cw"] = results["cw"] * INVPSI_TO_INVBAR  # 1/psi -> 1/bar
        # uw in cP — no conversion
        results["pb"] = results["pb"] * PSI_TO_BAR
        results["rsb"] = results["rsb"] * SCF_PER_STB_TO_SM3_PER_SM3

        # Undersaturated data (if PVTO)
        if pvto and results["usat"]:
            usat_p_out, usat_bo_out, usat_uo_out = results["usat"]
            for i in range(len(usat_p_out)):
                usat_p_out[i] = [p * PSI_TO_BAR for p in usat_p_out[i]]
                # usat_bo and usat_uo: Bo is dimensionless, uo is cP — no conversion
            results["usat"] = [usat_p_out, usat_bo_out, usat_uo_out]

    return results


def make_pvtw_table(
    pi: float,
    degf: float,
    wt: float = 0,
    ch4_sat: float = 0,
    pmin: float = 500,
    pmax: float = 10000,
    nrows: int = 20,
    export: bool = False,
    metric: bool = False,
) -> dict:
    """ Generates a PVTW (water PVT) table over a pressure range using brine_props (Spivey correlation).

        pi: Initial (reference) pressure (psia / barsa if metric)
        degf: Temperature (deg F / deg C if metric)
        wt: Salt wt% (0-100), default 0
        ch4_sat: Degree of methane saturation (0 - 1), default 0
        pmin: Minimum pressure for table (psia / barsa if metric), default 500
        pmax: Maximum pressure for table (psia / barsa if metric), default 10000
        nrows: Number of rows in table, default 20
        export: If True, writes PVTW.INC (ECLIPSE keyword) and pvtw_table.xlsx
        metric: If True, inputs/outputs use Eclipse METRIC units (barsa, degC, sm3/sm3, 1/bar).
                If False (default), inputs/outputs use FIELD units (psia, degF, scf/stb, 1/psi).

        Returns dict with keys:
            table: pandas DataFrame with columns Pressure, Bw, Density, Viscosity, Cw, Rsw
            pref: Reference pressure (psia / barsa)
            bw_ref: Bw at reference pressure (rb/stb / rm3/sm3)
            cw_ref: Compressibility at reference pressure (1/psi / 1/bar)
            visw_ref: Viscosity at reference pressure (cP)
            rsw_ref: Rsw at reference pressure (scf/stb / sm3/sm3)
            den_ref: Density (sg) at reference pressure
    """
    # --- Metric input conversion ---
    if metric:
        pi_orig_metric = pi
        degf_orig_metric = degf
        pi = pi * BAR_TO_PSI
        pmin = pmin * BAR_TO_PSI
        pmax = pmax * BAR_TO_PSI
        degf = degc_to_degf(degf)

    # Build pressure grid, ensuring pi is included (all in oilfield units internally)
    pressures = list(np.linspace(pmin, pmax, nrows))
    if pi not in pressures:
        pressures.append(pi)
    pressures = sorted(set(pressures))

    bws, dens, visws, cws_usat, cws_sat, rsws = [], [], [], [], [], []
    for p in pressures:
        bw, lden, visw, cw, rsw = brine.brine_props(p=p, degf=degf, wt=wt, ch4_sat=ch4_sat)
        bws.append(bw)
        dens.append(lden)
        visws.append(visw)
        cws_usat.append(cw[0])
        cws_sat.append(cw[1])
        rsws.append(rsw)

    df = pd.DataFrame()
    df["Pressure (psia)"] = pressures
    df["Bw (rb/stb)"] = bws
    df["Density (sg)"] = dens
    df["Viscosity (cP)"] = visws
    df["Cw_usat (1/psi)"] = cws_usat
    df["Cw_sat (1/psi)"] = cws_sat
    df["Rsw (scf/stb)"] = rsws

    # Reference properties at pi
    bw_ref, den_ref, visw_ref, cw_ref, rsw_ref = brine.brine_props(
        p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat
    )

    # --- Metric output conversion (done once; export reuses the converted frame) ---
    if metric:
        df["Pressure (psia)"] = df["Pressure (psia)"] * PSI_TO_BAR
        df.rename(columns={"Pressure (psia)": "Pressure (barsa)"}, inplace=True)
        # Bw is dimensionless FVF - no conversion
        df.rename(columns={"Bw (rb/stb)": "Bw (rm3/sm3)"}, inplace=True)
        # Density (sg) - no conversion
        # Viscosity (cP) - no conversion
        df["Cw_usat (1/psi)"] = df["Cw_usat (1/psi)"] * INVPSI_TO_INVBAR
        df.rename(columns={"Cw_usat (1/psi)": "Cw_usat (1/bar)"}, inplace=True)
        df["Cw_sat (1/psi)"] = df["Cw_sat (1/psi)"] * INVPSI_TO_INVBAR
        df.rename(columns={"Cw_sat (1/psi)": "Cw_sat (1/bar)"}, inplace=True)
        df["Rsw (scf/stb)"] = df["Rsw (scf/stb)"] * SCF_PER_STB_TO_SM3_PER_SM3
        df.rename(columns={"Rsw (scf/stb)": "Rsw (sm3/sm3)"}, inplace=True)

    if export:
        if metric:
            # Use metric values for ECLIPSE export
            pref_exp = pi * PSI_TO_BAR
            cw_ref_exp = cw_ref[0] * INVPSI_TO_INVBAR  # Undersaturated for PVTW keyword
            temp_str = f"{degf_orig_metric:.1f} deg C"
            unit_label = "METRIC"
        else:
            pref_exp = pi
            cw_ref_exp = cw_ref[0]  # Undersaturated for PVTW keyword
            temp_str = f"{degf:.1f} deg F"
            unit_label = "FIELD"

        # Write ECLIPSE PVTW keyword
        # PVTW format: Pref  Bw  Cw  Visw  Viscosibility
        # Viscosibility set to 0 (constant viscosity assumption)
        pvtw_line = f"  {pref_exp:.1f}  {bw_ref:.6f}  {cw_ref_exp:.6e}  {visw_ref:.4f}  0.0"
        fileout = f"-- Generated by pyResToolbox make_pvtw_table\n"
        fileout += f"-- Unit system: {unit_label}\n"
        fileout += f"-- Temperature: {temp_str}, Salt: {wt:.1f} wt%, CH4 sat: {ch4_sat:.2f}\n"
        fileout += f"PVTW\n{pvtw_line} /\n"
        with open("PVTW.INC", "w") as f:
            f.write(fileout)

        # Write full table to Excel (df already in output units)
        df.to_excel("pvtw_table.xlsx", index=False)

    if metric:
        return {
            "table": df,
            "pref": pi * PSI_TO_BAR,
            "bw_ref": bw_ref,
            "cw_ref": [c * INVPSI_TO_INVBAR for c in cw_ref],
            "visw_ref": visw_ref,
            "rsw_ref": rsw_ref * SCF_PER_STB_TO_SM3_PER_SM3,
            "den_ref": den_ref,
        }

    return {
        "table": df,
        "pref": pi,
        "bw_ref": bw_ref,
        "cw_ref": cw_ref,
        "visw_ref": visw_ref,
        "rsw_ref": rsw_ref,
        "den_ref": den_ref,
    }
