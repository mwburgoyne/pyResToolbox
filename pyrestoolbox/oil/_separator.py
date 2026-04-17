"""Separator gas and stock-tank gas utility functions."""

import numpy as np

from pyrestoolbox.constants import BAR_TO_PSI, degc_to_degf
from pyrestoolbox.shared_fns import validate_pe_inputs
from ._constants import _SGEVOL_THRESHOLD, _SGEVOL_HIGH, _SGEVOL_LOW


def sg_evolved_gas(
    p: float, degf: float, rsb: float, api: float, sg_sp: float
) -> float:
    """ Returns estimated specific gravity of gas evolved from oil insitu due to depressurization below Pb
        uses McCain & Hill Correlation (1995, SPE 30773)

        p: Pressure (psia)
        degf: Temperature (deg F)
        rsb: Oil solution GOR at Pb (scf/stb)
        api: Stock tank oil density (API)
        sg_sp: Specific gravity of separator gas (relative to air)
    """
    validate_pe_inputs(p=p, degf=degf, sg=sg_sp)

    if (
        p > _SGEVOL_THRESHOLD
    ):  # Two different sets from original 1995 paper (not reflected in Correlations book)
        a = _SGEVOL_HIGH
    else:
        a = _SGEVOL_LOW
    one_on_sgr = (
        a[1] / p
        + a[2] / p ** 2
        + a[3] * p
        + a[4] / degf ** 0.5
        + a[5] * degf
        + a[6] * rsb
        + a[7] * api
        + a[8] / sg_sp
        + a[9] * sg_sp ** 2
    )  # Eq 3.25
    return max(1 / one_on_sgr, sg_sp)


def sg_st_gas(
    psp: float, rsp: float, api: float, sg_sp: float, degf_sp: float
) -> float:
    """ Estimates specific gravity of gas evolving from stock tank
        from oil API and separator gas properties & conditions
        Returns sg_st (Stock Tank gas SG relative to air).
        Correlation reproduced from Valko McCain 2003 paper Eq 4-2

        psp: Separator pressure (psia)
        rsp: Separator GOR (separator scf / stb)
        api: Stock tank oil density (API)
        sg_sp: Separator gas specific gravity relative to air
        degf_sp: Separator temperature (deg f)
    """
    var = [np.log(psp), np.log(rsp), api, sg_sp, degf_sp]
    C = [
        [-17.275, -0.3354, 3.705, -155.52, 2.085],
        [7.9597, -0.3346, -0.4273, 629.61, -7.097e-2],
        [-1.1013, 0.1956, 1.818e-2, -957.38, 9.859e-4],
        [2.7735e-2, -3.4374e-2, -3.459e-4, 647.57, -6.312e-6],
        [3.2287e-3, 2.08e-3, 2.505e-6, -163.26, 1.4e-8],
    ]
    Zn = [sum([C[i][n] * var[n] ** i for i in range(5)]) for n in range(5)]
    Z = sum(Zn)
    sg_st = (
        1.219 + 0.198 * Z + 0.0845 * Z ** 2 + 0.03 * Z ** 3 + 0.003 * Z ** 4
    )
    return sg_st


def sgg_wt_avg(sg_sp: float, rsp: float, sg_st: float, rst: float) -> float:
    """ Calculates weighted average specific gravity of surface gas (sg_g)
        from separator and stock tank gas properties
        Returns sg_g (Weighted average surface gas SG relative to air).
        From McCain Correlations book, Eq 3.4

        sg_sp: Separator gas specific gravity relative to air
        rsp: Separator GOR (separator scf / stb)
        sg_st: Stock tank gas specific gravity relative to air
        rst: Stock tank producing gas-oil ratio (scf/stb)
    """
    sg_g = (sg_sp * rsp + sg_st * rst) / (rsp + rst)
    return sg_g


def oil_rs_st(psp: float, degf_sp: float, api: float, metric: bool = False) -> float:
    """ Estimates incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (scf/stb)

        Rsb = Rsp + Rst (Solution GOR at bubble point = Separator GOR + Stock Tank GOR).
        In absence of separator properties, a simple linear relationship with Rsp could be used instead;
          rs_st = 0.1618 * Separator GOR (Adapted from Eq 3-4 in Valko McCain 2003 paper)
        Correlation reproduced from Valko McCain 2003 paper Eq 3-2

        psp: Separator pressure (psia | barsa)
        degf_sp: Separator temperature (deg f | deg C)
        api: Stock tank oil density (API)
        metric: If True, input/output in Eclipse METRIC units. Defaults to False (FIELD)
    """
    if metric:
        psp = psp * BAR_TO_PSI
        degf_sp = degc_to_degf(degf_sp)

    validate_pe_inputs(p=psp, degf=degf_sp)

    var = [np.log(psp), np.log(degf_sp), api]
    C = [[-8.005, 1.224, -1.587], [2.7, -0.5, 0.0441], [-0.161, 0, -2.29e-5]]
    Zn = [sum([C[i][n] * var[n] ** i for i in range(3)]) for n in range(3)]
    Z = sum(Zn)
    return max(0, 3.955 + 0.83 * Z - 0.024 * Z ** 2 + 0.075 * Z ** 3)
