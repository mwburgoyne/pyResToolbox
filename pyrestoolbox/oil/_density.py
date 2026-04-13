"""Live oil density calculations."""

import warnings
import numpy as np

from pyrestoolbox.constants import (
    psc, tsc, BAR_TO_PSI, degc_to_degf,
    SM3_PER_SM3_TO_SCF_PER_STB, LBCUFT_TO_KGM3,
)
from pyrestoolbox.classes import deno_method
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.shared_fns import validate_pe_inputs
from pyrestoolbox._accelerator import RUST_AVAILABLE as _RUST_AVAILABLE
if _RUST_AVAILABLE:
    from pyrestoolbox import _native as _rust

from ._constants import (
    _API_NUMER, _API_DENOM, _COFB_C, _COFB_A0, _COFB_A1, _COFB_A2,
    _SWMH_RHOPO_INIT, _SWMH_RHOPO_RS, _SWMH_RHOA,
    _SWMH_MASS_NUMER_OIL, _SWMH_MASS_DENOM,
    _SWMH_SG_RHOA_A, _SWMH_SG_RHOA_B, _SWMH_SG_RHOA_C, _SWMH_SG_RHOA_D,
    _SWMH_DP_A, _SWMH_DP_B, _SWMH_DP_C, _SWMH_DP_D, _SWMH_DP_E, _SWMH_DP_F, _SWMH_DP_G,
    _SWMH_DT_A, _SWMH_DT_B, _SWMH_DT_C, _SWMH_DT_D, _SWMH_DT_E, _SWMH_DT_F, _SWMH_DT_G, _SWMH_DT_H,
)
from ._utils import check_sgs, oil_sg


def _cofb_mccain(api, sg_sp, pb, p, rsb, degf):
    """McCain Eq 3.13 cofb polynomial for undersaturated oil compressibility."""
    var = [
        np.log(api), np.log(sg_sp), np.log(pb),
        np.log(p / pb), np.log(rsb), np.log(degf),
    ]
    Zn = [sum([_COFB_C[i][n] * var[n] ** i for i in range(3)]) for n in range(6)]
    Zp = sum(Zn)
    ln_cofb = _COFB_A0 + _COFB_A1 * Zp + _COFB_A2 * Zp ** 2 - np.log(1e6)
    return np.exp(ln_cofb)


def oil_deno(
    p: float,
    degf: float,
    rs: float,
    rsb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    pb: float = 1e6,
    sg_o: float = 0,
    api: float = 0,
    denomethod: deno_method = deno_method.SWMH,
    metric: bool = False,
) -> float:
    """ Returns live oil density (lb/cuft | kg/m3) calculated with different correlations

        p: Pressure (psia | barsa)
        pb: Bubble point pressure (psia | barsa). Defaults to 1E6, and not used for densities below Pb. A valid value is required for density calculations above Pb
        degf: Reservoir Temperature (deg F | deg C)
        rs: Oil solution gas volume (scf/stb | sm3/sm3)
        rsb: Oil solution gas volume at bubble point pressure (scf/stb | sm3/sm3)
        sg_g: Weighted average specific gravity of surface gas (relative to air).
        sg_sp: Separator gas specific gravity (relative to air). If not known, an alternate nethod to estimate pseudo liquid density of surface gas will be used
        sg_o: Stock tank oil specific gravity (SG relative to water). If undefined will calculate from api
        api: Stock tank oil density (deg API). If undefined will calculate from sg_o. If both defined api value will prevail
        denomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   SWMH: Standing, Witte, McCain-Hill (1995) - Default
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3, kg/m3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        rs = rs * SM3_PER_SM3_TO_SCF_PER_STB
        rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        if pb != 1e6:
            pb = pb * BAR_TO_PSI

    validate_pe_inputs(p=p, degf=degf)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    denomethod = validate_methods(["denomethod"], [denomethod])

    if sg_g == 0 and sg_sp == 0:
        raise ValueError(
            "Must define at least one of sg_g and sg_sp for density calculation"
        )

    # Density at or below initial bubble point pressure
    def Deno_standing_white_mccainhill(
        p: float,
        degf: float,
        rs: float,
        rsb: float,
        sg_g: float,
        sg_sp: float,
        pb: float,
        sg_o: float,
        api: float,
    ) -> float:  # (1995), Eq 3.18a - 3.18g

        if sg_sp > 0:
            a = _SWMH_RHOA
            rho_po = max(_SWMH_RHOPO_INIT + _SWMH_RHOPO_RS * rs, 20)  # First estimate
            err = 1
            i = 0
            while err > 1e-8:
                i += 1
                rhoa = (
                    a[0]
                    + a[1] * sg_sp
                    + a[2] * sg_sp * rho_po
                    + a[3] * sg_sp * rho_po ** 2
                    + a[4] * rho_po
                    + a[5] * rho_po ** 2
                )  # Eq 3.18c
                new_rho_po = (rs * sg_sp + _SWMH_MASS_NUMER_OIL * sg_o) / (
                    _SWMH_MASS_DENOM + rs * sg_sp / rhoa
                )  # pseudoliquid density, Eq 3.18b. Note equation in original paper uses sg_sp rather than sg_g as in book.
                err = abs(rho_po - new_rho_po)
                rho_po = new_rho_po
                if i > 100:
                    warnings.warn(f"oil_deno: rho_po iteration did not converge after 100 iterations (err={err:.2e})")
                    break
        else:
            rhoa = _SWMH_SG_RHOA_A * (10 ** (_SWMH_SG_RHOA_B * api)) + (
                _SWMH_SG_RHOA_C + _SWMH_SG_RHOA_D * np.log10(api)
            ) * np.log10(
                sg_g
            )  # Eq 3.17e using sg_g. Apparent liquid density of surface gases
            rho_po = (rs * sg_g + _SWMH_MASS_NUMER_OIL * sg_o) / (
                _SWMH_MASS_DENOM + rs * sg_g / rhoa
            )  # pseudoliquid density, Eq 3.18b

        drho_p = (
            (_SWMH_DP_A + _SWMH_DP_B * 10 ** (_SWMH_DP_C * rho_po)) * p / 1000
            - _SWMH_DP_D * (_SWMH_DP_E + _SWMH_DP_F * 10 ** (_SWMH_DP_G * rho_po)) * (p / 1000) ** 2
        )  # Eq 3.19d

        rho_bs = rho_po + drho_p  # fake density used in calculations, Eq 3.19e
        drho_t = (
            (_SWMH_DT_A + _SWMH_DT_B * rho_bs ** _SWMH_DT_C) * (degf - tsc) ** _SWMH_DT_D
            - (_SWMH_DT_E + _SWMH_DT_F * 10 ** (_SWMH_DT_G * rho_bs))
            * (degf - tsc) ** _SWMH_DT_H
        )  # Eq 3.19f
        rho_or = rho_bs - drho_t  # Eq 3.19g

        return rho_or

    def Deno_p_gt_pb(
        p: float,
        degf: float,
        rs: float,
        rsb: float,
        sg_g: float,
        sg_sp: float,
        pb: float,
        sg_o: float,
        api: float,
    ) -> float:
        rhorb = Deno_standing_white_mccainhill(
            p=pb,
            degf=degf,
            rs=rs,
            rsb=rsb,
            sg_g=sg_g,
            sg_sp=sg_sp,
            pb=pb,
            sg_o=sg_o,
            api=api,
        )

        cofb_p = _cofb_mccain(api, sg_sp, pb, p, rsb, degf)
        return rhorb * np.exp(cofb_p * (p - pb))  # Eq 3.20

    fn_dic = {
        "SWMH": Deno_standing_white_mccainhill,
        "PGTPB": Deno_p_gt_pb,
    }  # Pressure greater than Pb

    if api == 0 and sg_o == 0:
        raise ValueError("Must supply either sg_o or api")

    if api == 0:  # Set api from sg_o
        api = _API_NUMER / sg_o - _API_DENOM
    else:  # overwrite sg_o with api value
        sg_o = oil_sg(api)

    if _RUST_AVAILABLE and denomethod.name == "SWMH":
        try:
            result = _rust.oil_deno_mccain_rust(p, degf, rs, rsb, sg_g, sg_sp, pb, sg_o, api)
            if metric:
                return result * LBCUFT_TO_KGM3
            return result
        except (ImportError, AttributeError):
            pass

    if (
        p > pb
    ):  # Use Eq 3.20, calculating oil density from density at Pb and compressibility factor
        result = fn_dic["PGTPB"](
            p=p,
            degf=degf,
            rs=rs,
            rsb=rsb,
            sg_g=sg_g,
            sg_sp=sg_sp,
            pb=pb,
            sg_o=sg_o,
            api=api,
        )
        if metric:
            return result * LBCUFT_TO_KGM3  # lb/cuft -> kg/m3
        return result

    result = fn_dic[denomethod.name](
        p=p,
        degf=degf,
        rs=rs,
        rsb=rsb,
        sg_g=sg_g,
        sg_sp=sg_sp,
        pb=pb,
        sg_o=sg_o,
        api=api,
    )
    if metric:
        return result * LBCUFT_TO_KGM3  # lb/cuft -> kg/m3
    return result
