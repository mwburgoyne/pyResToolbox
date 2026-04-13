"""Black oil table generation helpers."""

import numpy as np
import pandas as pd
from tabulate import tabulate

from pyrestoolbox.constants import psc, tsc, CUFTperBBL, WDEN
import pyrestoolbox.gas as gas
import pyrestoolbox.brine as brine

from ._utils import oil_sg
from ._correlations import oil_pbub, oil_rs_bub, oil_rs, oil_bo, oil_deno, oil_viso
from ._compressibility import oil_co
from ._harmonize import oil_harmonize


def _resolve_pb_rsb(pb, rsb, degf, api, sg_sp, sg_g, pvto, pmax,
                     rsmethod, pbmethod):
    """Internal helper: resolves Pb/Rsb using oil_harmonize plus PVTO extension.

    Returns (pb, rsb, rsb_frac, rsb_max, pb_i, rsb_i).
    """
    pb_i = pb
    rsb_i = rsb

    pb, rsb, rsb_frac, _vis_frac = oil_harmonize(
        pb=pb, rsb=rsb, degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
        rsmethod=rsmethod, pbmethod=pbmethod,
    )

    rsb_max = rsb

    if pvto and pmax > pb:
        rsb_max = oil_rs_bub(
            degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
            pb=pmax, rsmethod=rsmethod,
        )
        rs_at_pbi = oil_rs(
            api=api, degf=degf, sg_sp=sg_sp, p=pb, pb=pmax,
            rsb=rsb_max, rsmethod=rsmethod, pbmethod=pbmethod,
        )
        err = rs_at_pbi - rsb
        rsb_frac_new = rsb_frac
        i = 0
        while abs(err) > 0.0001:
            rsb_frac_new = rsb / rs_at_pbi * rsb_frac
            rs_at_pbi = oil_rs(
                api=api, degf=degf, sg_sp=sg_sp, p=pb, pb=pmax,
                rsb=rsb_max * rsb_frac_new,
                rsmethod=rsmethod, pbmethod=pbmethod,
            )
            rsb_frac = rsb_frac_new
            err = rs_at_pbi - rsb
            i += 1
            if i > 100:
                raise RuntimeError(
                    "Could not solve Pb & Rsb for these combination of inputs"
                )
        rsb_frac = rsb_frac_new

    return pb, rsb, rsb_frac, rsb_max, pb_i, rsb_i


def _build_bot_tables(pressures, pb, rsb, rsb_frac, rsb_max, sg_o, sg_g, sg_sp,
                      api, degf, pvto, wt, ch4_sat,
                      zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod,
                      vis_frac=1.0):
    """Compute all PVT properties over the pressure array.

    Returns (rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
             usat_p, usat_bo, usat_uo) where usat_* are empty lists if not pvto.
    """
    if pvto:
        pb = max(pb, max(pressures))
        rsb = rsb_max * rsb_frac

    co, cg, rss, bos, uos, gfvf, visg, gz, bws, visws, denos = [
        [] for _ in range(11)
    ]

    for p in pressures:
        if p > pb:
            rss.append(rsb)
        else:
            rss.append(
                oil_rs(
                    api=api, degf=degf, sg_sp=sg_sp, p=p, pb=pb,
                    rsb=rsb / rsb_frac, rsmethod=rsmethod, pbmethod=pbmethod,
                ) * rsb_frac
            )

        bos.append(
            oil_bo(
                p=p, pb=pb, degf=degf, rs=rss[-1], rsb=rsb,
                sg_g=sg_g, sg_sp=sg_sp, sg_o=sg_o,
                denomethod=denomethod, bomethod=bomethod,
            )
        )
        denos.append(
            oil_deno(
                p=p, degf=degf, rs=rss[-1], rsb=rsb,
                sg_g=sg_g, sg_sp=sg_sp, pb=pb, sg_o=sg_o, api=api,
            )
        )
        uos.append(oil_viso(p=p, api=api, degf=degf, pb=pb, rs=rss[-1]) * vis_frac)
        co.append(
            oil_co(
                p=p, api=api, sg_sp=sg_sp, sg_g=sg_g, degf=degf,
                pb=pb, rsb=rsb, undersaturated_only=True,
                zmethod=zmethod, rsmethod=rsmethod,
                cmethod=cmethod, denomethod=denomethod, bomethod=bomethod,
            )
        )

        gfvf.append(
            gas.gas_bg(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
            * 1000 / CUFTperBBL
        )
        gz.append(
            gas.gas_z(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        visg.append(
            gas.gas_ug(p=p, sg=sg_g, degf=degf, zmethod=zmethod, cmethod=cmethod)
        )
        cg.append(gas.gas_cg(p=p, sg=sg_g, degf=degf, cmethod=cmethod))
        bw, _lden, visw, _cw, _rsw = brine.brine_props(
            p=p, degf=degf, wt=wt, ch4_sat=ch4_sat
        )
        bws.append(bw)
        visws.append(visw)

    # Undersaturated extension for PVTO
    usat_p, usat_bo, usat_uo = [], [], []
    if pvto:
        for i, p in enumerate(pressures):
            if i == 0:
                continue
            try:
                usat_p.append(pressures[i:])
                usat_bo.append(
                    [
                        oil_bo(
                            p=pusat, pb=p, degf=degf, rs=rss[i], rsb=rss[i],
                            sg_g=sg_g, sg_sp=sg_sp, sg_o=sg_o,
                            denomethod=denomethod, bomethod=bomethod,
                        )
                        for pusat in usat_p[-1]
                    ]
                )
                usat_uo.append(
                    [
                        oil_viso(p=pusat, api=api, degf=degf, pb=p, rs=rss[i]) * vis_frac
                        for pusat in usat_p[-1]
                    ]
                )
            except (ValueError, IndexError, ZeroDivisionError):
                pass

    return (rss, bos, denos, uos, co, gz, gfvf, cg, visg, bws, visws,
            usat_p, usat_bo, usat_uo)


def _format_bot_results(pressures, rss, bos, denos, uos, co, gz, gfvf, cg,
                        visg, bws, visws, usat_p, usat_bo, usat_uo,
                        sg_o, sg_g, pi, degf, wt, ch4_sat, pb_i, rsb_i,
                        rsb_frac, pvto, export, zmethod, cmethod,
                        vis_frac=1.0):
    """Assemble DataFrame, optionally export Eclipse files, and return results dict."""
    st_deno = sg_o * WDEN
    st_deng = gas.gas_den(
        p=psc, sg=sg_g, degf=tsc, zmethod=zmethod, cmethod=cmethod
    )
    bw, lden, visw, cw_list, _rsw = brine.brine_props(
        p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat
    )
    res_denw = lden * WDEN
    res_cw = cw_list[0]  # Undersaturated compressibility for BOT

    df = pd.DataFrame()
    df["Pressure (psia)"] = pressures
    df["Rs (mscf/stb)"] = rss
    df["Rs (mscf/stb)"] = df["Rs (mscf/stb)"] / 1000
    df["Bo (rb/stb)"] = bos
    df["Deno (lb/cuft)"] = denos
    df["uo (cP)"] = uos
    df["Co (1/psi)"] = co
    df["Gas Z (v/v)"] = gz
    df["Bg (rb/mscf"] = gfvf
    df["Cg (1/psi)"] = cg
    df["ug (cP)"] = visg
    df["Bw (rb/stb)"] = bws
    df["uw (cP)"] = visws

    if export:
        df.to_excel("bot.xlsx", index=False, engine="openpyxl")
        pvdg = df[["Pressure (psia)", "Bg (rb/mscf", "ug (cP)"]]
        pvdg = pvdg.set_index("Pressure (psia)")
        headers = ["-- P (psia)", "Bg (rb/mscf", "ug (cP)"]
        fileout = "PVDG\n" + tabulate(pvdg, headers) + "\n/"
        with open("PVDG.INC", "w") as text_file:
            text_file.write(fileout)
        pvdo = df[["Pressure (psia)", "Bo (rb/stb)", "uo (cP)"]]
        pvdo = pvdo.set_index("Pressure (psia)")
        headers = ["-- P (psia)", "Bo (rb/stb)", "uo (cP)"]
        fileout = "PVDO\n" + tabulate(pvdo, headers) + "\n/"
        with open("PVDO.INC", "w") as text_file:
            text_file.write(fileout)

        if pvto:
            pvto_out = "PVTO\n"
            headers = [
                "-- Rs (mscf/stb)",
                "P (psia)",
                "Bo (rb/stb)",
                "uo (cP)",
                "",
            ]
            table = []
            for r, row in df.iterrows():
                table.append(
                    [
                        row["Rs (mscf/stb)"],
                        row["Pressure (psia)"],
                        row["Bo (rb/stb)"],
                        row["uo (cP)"],
                        "/",
                    ]
                )
                try:
                    if r > 0:
                        for e, entry in enumerate(usat_p[r - 1]):
                            if e == 0:
                                continue
                            table.append(
                                [
                                    " ",
                                    entry,
                                    usat_bo[r - 1][e],
                                    usat_uo[r - 1][e],
                                    " ",
                                ]
                            )
                except (IndexError, KeyError):
                    pass
            pvto_out += tabulate(table, headers)
            pvto_out += "\n/"
            with open("PVTO.INC", "w") as text_file:
                text_file.write(pvto_out)

    results = {
        "bot": df,
        "deno": st_deno,
        "deng": st_deng,
        "denw": res_denw,
        "cw": res_cw,
        "uw": visw,
        "pb": pb_i,
        "rsb": rsb_i,
        "rsb_scale": rsb_frac,
        "vis_frac": vis_frac,
        "usat": [],
    }
    if pvto:
        results["usat"] = [usat_p, usat_bo, usat_uo]

    return results


def make_bot_og(*args, **kwargs):
    """Deprecated: Use simtools.make_bot_og() instead. This wrapper remains for backward compatibility."""
    from pyrestoolbox.simtools import make_bot_og as _make_bot_og
    return _make_bot_og(*args, **kwargs)
