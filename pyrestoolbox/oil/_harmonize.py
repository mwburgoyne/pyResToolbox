"""Oil PVT harmonization: resolve consistent Pb, Rsb, rsb_frac, vis_frac."""

from typing import Tuple

from pyrestoolbox.constants import (
    psc, BAR_TO_PSI, PSI_TO_BAR, degc_to_degf,
    SM3_PER_SM3_TO_SCF_PER_STB, SCF_PER_STB_TO_SM3_PER_SM3,
)
from pyrestoolbox.classes import pb_method, rs_method
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.shared_fns import validate_pe_inputs

from ._utils import check_sgs
from ._correlations import oil_pbub, oil_rs_bub, oil_rs, oil_viso


def oil_harmonize(
    pb: float = 0,
    rsb: float = 0,
    degf: float = 0,
    api: float = 0,
    sg_sp: float = 0,
    sg_g: float = 0,
    uo_target: float = 0,
    p_uo: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
    pbmethod: pb_method = pb_method.VELAR,
    metric: bool = False,
) -> Tuple:
    """Resolves consistent Pb, Rsb, rsb_frac, and vis_frac from user inputs.

    Given one or both of Pb and Rsb, returns values that are mutually consistent
    with the selected oil PVT correlations. Optionally computes a viscosity scaling
    factor (vis_frac) to match a known viscosity measurement.

    - If only pb is specified (rsb=0): calculates rsb from pb
    - If only rsb is specified (pb=0): calculates pb from rsb
    - If both are specified: finds rsb_frac scaling factor that honors both values
    - If uo_target and p_uo are specified: computes vis_frac = uo_target / uo_corr

    pb: Bubble point pressure (psia | barsa). Default 0 (unknown)
    rsb: Solution GOR at Pb (scf/stb | sm3/sm3). Default 0 (unknown)
    degf: Reservoir temperature (deg F | deg C)
    api: Stock tank oil density (deg API)
    sg_sp: Separator gas specific gravity
    sg_g: Weighted average surface gas specific gravity
    uo_target: Target oil viscosity (cP) at pressure p_uo. Default 0 (no viscosity tuning)
    p_uo: Pressure at which uo_target was measured (psia | barsa). Required if uo_target > 0
    rsmethod: Rs calculation method. Default VELAR
    pbmethod: Pb calculation method. Default VELAR
    metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)

    Returns tuple of (pb, rsb, rsb_frac, vis_frac) where:
      - pb: Bubble point pressure (psia | barsa)
      - rsb: Solution GOR at bubble point (scf/stb | sm3/sm3)
      - rsb_frac: Scaling factor applied to Rs correlation to honor user Pb and Rsb
                  (1.0 if only one was specified)
      - vis_frac: Multiplicative viscosity scaling factor (1.0 if no target specified)
    """
    if metric:
        degf = degc_to_degf(degf)
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        if p_uo > 0:
            p_uo = p_uo * BAR_TO_PSI

    validate_pe_inputs(degf=degf)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    rsmethod, pbmethod = validate_methods(
        ["rsmethod", "pbmethod"], [rsmethod, pbmethod]
    )

    rsb_i = rsb
    pb_i = pb
    rsb_frac = 1.0

    # Calculate rsb from pb
    if rsb_i <= 0 and pb_i > 0:
        rsb = oil_rs_bub(
            degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
            pb=pb, rsmethod=rsmethod,
        )

    # Calculate pb from rsb
    if pb_i <= 0 and rsb_i > 0:
        pb = oil_pbub(
            degf=degf, api=api, sg_sp=sg_sp, rsb=rsb, pbmethod=pbmethod
        )

    # Both defined by user — find rsb_frac that honors both
    if pb_i > 0 and rsb_i > 0:
        pbcalc = oil_pbub(
            degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
            rsb=rsb, pbmethod=pbmethod,
        )
        err = 100
        rsb_old = rsb
        i = 0
        while err > 0.0001:
            rsbnew = pb / pbcalc * rsb_old
            pbcalc = oil_pbub(
                degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
                rsb=rsbnew, pbmethod=pbmethod,
            )
            rsb_old = rsbnew
            err = abs(pb - pbcalc)
            i += 1
            if i > 100:
                raise RuntimeError(
                    "Could not solve Pb & Rsb for these combination of inputs"
                )
        rsb_frac = rsb_i / rsbnew

    # Compute vis_frac if target viscosity specified
    vis_frac = 1.0
    if uo_target > 0 and p_uo > 0:
        rs_at_p = oil_rs(
            api=api, degf=degf, sg_sp=sg_sp, p=p_uo,
            pb=pb, rsb=rsb / rsb_frac, rsmethod=rsmethod, pbmethod=pbmethod,
        ) * rsb_frac
        uo_corr = oil_viso(p=p_uo, api=api, degf=degf, pb=pb, rs=rs_at_p)
        vis_frac = uo_target / uo_corr

    if metric:
        return pb * PSI_TO_BAR, rsb * SCF_PER_STB_TO_SM3_PER_SM3, rsb_frac, vis_frac
    return pb, rsb, rsb_frac, vis_frac


def oil_harmonize_pb_rsb(*args, **kwargs) -> Tuple:
    """Deprecated: Use oil_harmonize() instead.

    Returns (pb, rsb, rsb_frac) -- the original 3-tuple without vis_frac.
    """
    pb, rsb, rsb_frac, _vis_frac = oil_harmonize(*args, **kwargs)
    return pb, rsb, rsb_frac
