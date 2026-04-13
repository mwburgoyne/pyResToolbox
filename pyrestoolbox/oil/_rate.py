"""Radial and linear oil flow rate calculations."""

import numpy as np
import numpy.typing as npt

from pyrestoolbox.constants import (
    BAR_TO_PSI, degc_to_degf, M_TO_FT, SQM_TO_SQFT, STB_TO_SM3,
)
from pyrestoolbox.shared_fns import validate_pe_inputs

from ._constants import (
    _DARCY_RADIAL, _DARCY_SKIN_OFFSET,
    _VOGEL_AOF_DENOM, _VOGEL_LIN, _VOGEL_QUAD,
)
from ._correlations import oil_bo, oil_viso


def oil_rate_radial(
    k: npt.ArrayLike,
    h: npt.ArrayLike,
    pr: npt.ArrayLike,
    pwf: npt.ArrayLike,
    r_w: float,
    r_ext: float,
    uo: float = 0,
    bo: float = 0,
    S: float = 0,
    vogel: bool = False,
    pb: float = 0,
    oil_pvt = None,
    degf: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns liquid rate for radial flow (stb/day | sm3/day) using Darcy pseudo steady state equation
        k: Effective Permeability to flow (mD)
        h: Net flow height (ft | m)
        Pr: Reservoir pressure (psia | barsa)
        pwf: BHFP (psia | barsa)
        r_w: Wellbore Radius (ft | m)
        r_ext: External Reservoir Radius (ft | m)
        uo: Liquid viscosity (cP). Not required if oil_pvt is provided
        bo: Liquid Formation Volume Factor (rb/stb | rm3/sm3). Not required if oil_pvt is provided
        S: Wellbore Skin (Dimensionless). Defaults to zero if not specified
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
        pb: Bubble point pressure (psia | barsa). Defaults to zero if not specified. Not used unless Vogel option is invoked
        oil_pvt: OilPVT object. If provided, uo and bo are calculated from the PVT object at reservoir pressure
        degf: Reservoir temperature (deg F | deg C). Required when oil_pvt is provided
        metric: If True, input/output in Eclipse METRIC units (barsa, m, sm3/d). Defaults to False (FIELD)
    """
    if metric:
        h = np.asarray(h) * M_TO_FT if not isinstance(h, (int, float)) else h * M_TO_FT
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        r_w = r_w * M_TO_FT
        r_ext = r_ext * M_TO_FT
        if pb > 0:
            pb = pb * BAR_TO_PSI

    validate_pe_inputs(p=pr)
    validate_pe_inputs(p=pwf)
    if degf > 0:
        validate_pe_inputs(degf=degf)

    k, h, pr, pwf = (
        np.asarray(k),
        np.asarray(h),
        np.asarray(pr),
        np.asarray(pwf),
    )

    if oil_pvt is not None:
        if degf <= 0:
            raise ValueError("degf required when using oil_pvt")
        degf_f = degc_to_degf(degf) if metric else degf
        p_rep = float(np.asarray(pr).ravel()[0])
        rs_rep = oil_pvt._rs_field(p_rep, degf_f)
        uo = oil_viso(p=p_rep, api=oil_pvt.api, degf=degf_f, pb=oil_pvt.pb, rs=rs_rep) * oil_pvt.vis_frac
        bo = oil_bo(p=p_rep, pb=oil_pvt.pb, degf=degf_f, rs=rs_rep, rsb=oil_pvt.rsb,
                    sg_o=oil_pvt.sg_o, sg_g=oil_pvt.sg_g, sg_sp=oil_pvt.sg_sp,
                    bomethod=oil_pvt.bomethod)
        pb = oil_pvt.pb
        vogel = True
    elif uo <= 0 or bo <= 0:
        raise ValueError("Either oil_pvt or both uo and bo must be specified")

    if pwf.size > 1:
        pb = np.array([max(pb, pwf[i]) for i in range(pwf.size)])

    if r_w <= 0:
        raise ValueError("Wellbore radius r_w must be positive")
    if r_ext <= r_w:
        raise ValueError("External radius r_ext must be greater than wellbore radius r_w")

    J = (
        _DARCY_RADIAL * k * h / (uo * bo * (np.log(r_ext / r_w) + S - _DARCY_SKIN_OFFSET))
    )  # Productivity index
    if not vogel:
        qoil = J * (pr - pwf)
    else:
        if np.any(pb <= 0):
            raise ValueError("Bubble point pressure pb must be positive when vogel=True")
        qsat_max = J * pb / _VOGEL_AOF_DENOM
        qusat = J * (pr - pb)
        qoil = (
            qsat_max * (1 - _VOGEL_LIN * (pwf / pb) - _VOGEL_QUAD * (pwf / pb) ** 2) + qusat
        )
    if metric:
        return qoil * STB_TO_SM3  # stb/d -> sm3/d
    return qoil

def oil_rate_linear(
    k: npt.ArrayLike,
    pr: npt.ArrayLike,
    pwf: npt.ArrayLike,
    area: npt.ArrayLike,
    length: float,
    uo: float = 0,
    bo: float = 0,
    vogel: bool = False,
    pb: float = 0,
    oil_pvt = None,
    degf: float = 0,
    metric: bool = False,
) -> np.ndarray:
    """ Returns liquid rate for linear flow (stb/day | sm3/day) using Darcy steady state equation
        k: Permeability (mD)
        Pr: Reservoir pressure (psia | barsa)
        pwf: BHFP (psia | barsa)
        area: Net cross sectional area perpendicular to direction of flow (ft2 | m2).
        length: Length over which flow takes place (ft | m)
        uo: Liquid viscosity (cP). Not required if oil_pvt is provided
        bo: Liquid Formation Volume Factor (rb/stb | rm3/sm3). Not required if oil_pvt is provided
        vogel: (True / False). Invokes the Vogel model that reduces inflow below bubble point pressure. Defaults to False if undefined
        pb: Bubble point pressure (psia | barsa). Defaults to zero if not specified. Not used unless Vogel option is invoked
        oil_pvt: OilPVT object. If provided, uo and bo are calculated from the PVT object at reservoir pressure
        degf: Reservoir temperature (deg F | deg C). Required when oil_pvt is provided
        metric: If True, input/output in Eclipse METRIC units (barsa, m, m2, sm3/d). Defaults to False (FIELD)
    """
    if metric:
        pr = np.asarray(pr) * BAR_TO_PSI if not isinstance(pr, (int, float)) else pr * BAR_TO_PSI
        pwf = np.asarray(pwf) * BAR_TO_PSI if not isinstance(pwf, (int, float)) else pwf * BAR_TO_PSI
        area = np.asarray(area) * SQM_TO_SQFT if not isinstance(area, (int, float)) else area * SQM_TO_SQFT
        length = length * M_TO_FT
        if pb > 0:
            pb = pb * BAR_TO_PSI

    validate_pe_inputs(p=pr)
    validate_pe_inputs(p=pwf)
    if degf > 0:
        validate_pe_inputs(degf=degf)

    k, area, pr, pwf = (
        np.asarray(k),
        np.asarray(area),
        np.asarray(pr),
        np.asarray(pwf),
    )

    if oil_pvt is not None:
        if degf <= 0:
            raise ValueError("degf required when using oil_pvt")
        degf_f = degc_to_degf(degf) if metric else degf
        p_rep = float(np.asarray(pr).ravel()[0])
        rs_rep = oil_pvt._rs_field(p_rep, degf_f)
        uo = oil_viso(p=p_rep, api=oil_pvt.api, degf=degf_f, pb=oil_pvt.pb, rs=rs_rep) * oil_pvt.vis_frac
        bo = oil_bo(p=p_rep, pb=oil_pvt.pb, degf=degf_f, rs=rs_rep, rsb=oil_pvt.rsb,
                    sg_o=oil_pvt.sg_o, sg_g=oil_pvt.sg_g, sg_sp=oil_pvt.sg_sp,
                    bomethod=oil_pvt.bomethod)
        pb = oil_pvt.pb
        vogel = True
    elif uo <= 0 or bo <= 0:
        raise ValueError("Either oil_pvt or both uo and bo must be specified")

    if length <= 0:
        raise ValueError("Flow length must be positive")

    J = (
        _DARCY_RADIAL * k * area / (2 * np.pi * uo * bo * length)
    )  # Productivity index (linear Darcy: 0.00708/(2*pi) = 0.001127)
    if not vogel:
        qoil = J * (pr - pwf)
    else:
        if np.any(pb <= 0):
            raise ValueError("Bubble point pressure pb must be positive when vogel=True")
        qsat_max = J * pb / _VOGEL_AOF_DENOM
        qusat = J * (pr - pb)
        qoil = (
            qsat_max * (1 - _VOGEL_LIN * (pwf / pb) - _VOGEL_QUAD * (pwf / pb) ** 2) + qusat
        )
    if metric:
        return qoil * STB_TO_SM3  # stb/d -> sm3/d
    return qoil
