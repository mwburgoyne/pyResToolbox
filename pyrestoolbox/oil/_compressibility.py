"""Oil compressibility and total two-phase FVF."""

from pyrestoolbox.constants import (
    psc, CUFTperBBL, BAR_TO_PSI, degc_to_degf,
    SM3_PER_SM3_TO_SCF_PER_STB,
    INVPSI_TO_INVBAR,
)
from pyrestoolbox.classes import (
    z_method, c_method, pb_method, rs_method, bo_method, deno_method, co_method,
)
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.shared_fns import validate_pe_inputs
import pyrestoolbox.gas as gas

from ._utils import check_sgs, oil_sg
from ._density import _cofb_mccain
from ._correlations import oil_pbub, oil_rs_bub, oil_rs, oil_bo


def _perrine_co_sat(p, api, degf, sg_sp, sg_g, pb, rsb, zmethod, cmethod,
                    rsmethod, pbmethod, bomethod, denomethod):
    """Compute Perrine's saturated compressibility: co_sat = -(1/Bo)*dBo/dp + (Bg/Bo)*dRs/dp."""
    sg_o = oil_sg(api)
    dp = max(0.5, p * 0.001)
    p_hi = p + dp
    p_lo = max(p - dp, psc)
    span = p_hi - p_lo
    if span < 1e-10:
        span = dp

    rs_at_p = oil_rs(api=api, degf=degf, sg_sp=sg_sp, p=p, pb=pb, rsb=rsb,
                     rsmethod=rsmethod, pbmethod=pbmethod)
    bo_at_p = oil_bo(p=p, pb=pb, degf=degf, rs=rs_at_p, rsb=rsb, sg_sp=sg_sp,
                     sg_g=sg_g, sg_o=sg_o, bomethod=bomethod, denomethod=denomethod)
    rs_hi = oil_rs(api=api, degf=degf, sg_sp=sg_sp, p=p_hi, pb=pb, rsb=rsb,
                   rsmethod=rsmethod, pbmethod=pbmethod)
    rs_lo = oil_rs(api=api, degf=degf, sg_sp=sg_sp, p=p_lo, pb=pb, rsb=rsb,
                   rsmethod=rsmethod, pbmethod=pbmethod)
    bo_hi = oil_bo(p=p_hi, pb=pb, degf=degf, rs=rs_hi, rsb=rsb, sg_sp=sg_sp,
                   sg_g=sg_g, sg_o=sg_o, bomethod=bomethod, denomethod=denomethod)
    bo_lo = oil_bo(p=p_lo, pb=pb, degf=degf, rs=rs_lo, rsb=rsb, sg_sp=sg_sp,
                   sg_g=sg_g, sg_o=sg_o, bomethod=bomethod, denomethod=denomethod)
    dBodp = (bo_hi - bo_lo) / span
    dRsdp = (rs_hi - rs_lo) / span
    bg_at_p = gas.gas_bg(p, sg_sp, degf, zmethod=zmethod, cmethod=cmethod) / CUFTperBBL
    return -1.0 / bo_at_p * dBodp + bg_at_p / bo_at_p * dRsdp


def oil_co(
    p: float,
    api: float,
    degf: float,
    sg_sp: float = 0,
    sg_g: float = 0,
    pb: float = 0,
    rsb: float = 0,
    pi: float = 0,
    co_sat: bool = False,
    undersaturated_only: bool = False,
    comethod: co_method = co_method.EXPLT,
    zmethod: z_method = z_method.DAK,
    rsmethod: rs_method = rs_method.VELAR,
    cmethod: c_method = c_method.PMC,
    denomethod: deno_method = deno_method.SWMH,
    bomethod: bo_method = bo_method.MCAIN,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
) -> float:
    """ Returns oil compressibility (1/psi | 1/bar).

        By default (co_sat=False) returns undersaturated compressibility calculated with
        Co = -1/Bo * dBo/dp at constant Rs, using numerically derived values. Rs is held at the
        equilibrium value for the specified pressure (rsb above Pb, correlation value below Pb).

        When co_sat=True, returns [co_usat, co_sat] list. The saturated compressibility uses
        Perrine's definition: co_sat = -(1/Bo)*dBo/dp + (Bg/Bo)*dRs/dp, where both Bo and Rs
        vary with pressure. Above Pb, co_sat equals co_usat (no gas evolution).

        When undersaturated_only=True, uses the analytical cofb correlation (McCain Eq 3.13)
        at all pressures, giving a smooth curve with no discontinuity at Pb. This is the
        compressibility at constant composition (rsb) regardless of pressure relative to Pb.

        p: Reservoir pressure (psia | barsa)
        api: Stock tank oil density (deg API)
        sg_sp: Separator Gas specific Gravity (relative to air). If not defined, will use sg_g instead
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not defined, will use sg_sp instead
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point pressure (psia | barsa). If not provided, will attempt to calculate with Valko-McCain Pb Correlation
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb | sm3/sm3)
        co_sat: If True, return [co_usat, co_sat] list. Default False (returns float)
        undersaturated_only: If True, use analytical cofb correlation at all pressures (no Pb discontinuity). Default False
        comethod: A string or co_method Enum class that specifies calculation method for compressibility (currently only one option)
        zmethod: A string or z_method Enum class that specifies calculation method for gas Z-factor
        rsmethod: A string or rs_method Enum class that specifies calculation method for GOR
        cmethod: A string or c_method Enum class that specifies calculation method for gas critical properties
        denomethod: A string or deno_method Enum class that specifies calculation method for live oil density
        bomethod: A string or bo_method Enum class that specifies calculation method for oil FVF
        pbmethod: A string or pb_method Enum class that specifies calculation method for bubble point pressure
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3, 1/bar). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        if pi > 0:
            pi = pi * BAR_TO_PSI

    validate_pe_inputs(p=p, degf=degf)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    (
        zmethod,
        rsmethod,
        cmethod,
        denomethod,
        bomethod,
        pbmethod,
        comethod,
    ) = validate_methods(
        [
            "zmethod",
            "rsmethod",
            "cmethod",
            "denomethod",
            "bomethod",
            "pbmethod",
            "comethod",
        ],
        [zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod, comethod],
    )

    if pb <= 0:  # Calculate Pb
        pb = oil_pbub(
            api=api, degf=degf, rsb=rsb, sg_sp=sg_sp, pbmethod=pbmethod
        )
    if rsb <= 0:  # Calculate rsb
        rsb = oil_rs_bub(
            api=api,
            degf=degf,
            pb=pb,
            sg_sp=sg_sp,
            rsmethod=rsmethod,
        )

    # Analytical cofb compressibility (McCain Eq 3.13) — smooth across Pb.
    # When undersaturated_only=True, uses the cofb polynomial at all pressures
    # giving a continuous curve with no discontinuity at Pb. This is the
    # compressibility at constant composition (rsb), suitable for BOT tables.
    if undersaturated_only:
        co_usat = _cofb_mccain(api, sg_sp, pb, p, rsb, degf)

        if not co_sat:
            if metric:
                return co_usat * INVPSI_TO_INVBAR
            return co_usat
        if p >= pb:
            if metric:
                return [co_usat * INVPSI_TO_INVBAR, co_usat * INVPSI_TO_INVBAR]
            return [co_usat, co_usat]
        co_sat_val = _perrine_co_sat(p, api, degf, sg_sp, sg_g, pb, rsb, zmethod, cmethod,
                                     rsmethod, pbmethod, bomethod, denomethod)
        if metric:
            return [co_usat * INVPSI_TO_INVBAR, co_sat_val * INVPSI_TO_INVBAR]
        return [co_usat, co_sat_val]

    def Co_explicit(
        p=p,
        api=api,
        sg_sp=sg_sp,
        sg_g=sg_g,
        degf=degf,
        pb=pb,
        rsb=rsb,
        zmethod=zmethod,
        rsmethod=rsmethod,
        cmethod=cmethod,
        denomethod=denomethod,
        bomethod=bomethod,
    ):  # Explicit - Calculate with numerical derivatives
        # Undersaturated compressibility: co = -1/bo * dBo/dp at constant Rs
        # Rs is held constant at the equilibrium value for pressure p,
        # so the derivative captures only liquid-phase compression
        if p > pb:
            rs_fixed = rsb
        else:
            rs_fixed = oil_rs(
                api=api,
                degf=degf,
                sg_sp=sg_sp,
                p=p,
                pb=pb,
                rsb=rsb,
                rsmethod=rsmethod,
                pbmethod=pbmethod,
            )

        def calc_bo_at_p(p_eval):
            sg_o = oil_sg(api)
            return oil_bo(
                p=p_eval,
                pb=pb,
                degf=degf,
                rs=rs_fixed,
                rsb=rsb,
                sg_sp=sg_sp,
                sg_g=sg_g,
                sg_o=sg_o,
                bomethod=bomethod,
                denomethod=denomethod,
            )

        dp = max(0.5, p * 0.001)  # Relative step size for numerical derivative
        p_hi = p + dp
        p_lo = p - dp

        # Clamp to avoid negative pressures
        p_lo = max(p_lo, psc)

        span = p_hi - p_lo
        if span < 1e-10:
            span = dp  # fallback for degenerate cases

        dbodp = (calc_bo_at_p(p_hi) - calc_bo_at_p(p_lo)) / span
        bo = calc_bo_at_p(p)
        return -1 / bo * dbodp

    fn_dic = {"EXPLT": Co_explicit}

    co_usat = fn_dic[comethod.name](
        p=p,
        api=api,
        sg_sp=sg_sp,
        sg_g=sg_g,
        degf=degf,
        pb=pb,
        rsb=rsb,
        zmethod=zmethod,
        rsmethod=rsmethod,
        cmethod=cmethod,
        denomethod=denomethod,
        bomethod=bomethod,
    )

    if not co_sat:
        if metric:
            return co_usat * INVPSI_TO_INVBAR  # 1/psi -> 1/bar
        return co_usat

    if p >= pb:
        co_sat_val = co_usat
    else:
        co_sat_val = _perrine_co_sat(p, api, degf, sg_sp, sg_g, pb, rsb, zmethod, cmethod,
                                     rsmethod, pbmethod, bomethod, denomethod)

    if metric:
        return [co_usat * INVPSI_TO_INVBAR, co_sat_val * INVPSI_TO_INVBAR]
    return [co_usat, co_sat_val]

def oil_bt(
    p: float,
    api: float,
    degf: float,
    sg_sp: float = 0,
    sg_g: float = 0,
    pb: float = 0,
    rsb: float = 0,
    rsi: float = 0,
    zmethod: z_method = z_method.DAK,
    rsmethod: rs_method = rs_method.VELAR,
    cmethod: c_method = c_method.PMC,
    denomethod: deno_method = deno_method.SWMH,
    bomethod: bo_method = bo_method.MCAIN,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
) -> float:
    """ Returns total two-phase oil formation volume factor Bt (rb/stb | rm3/sm3).

        Bt = Bo + (Rsi - Rs) * Bg

        Above Pb, Rs = Rsi so Bt = Bo. Below Pb, Bt accounts for the reservoir
        volume of both the liquid oil and the gas that has evolved from it relative
        to the original solution GOR.

        p: Reservoir pressure (psia | barsa)
        api: Stock tank oil density (deg API)
        sg_sp: Separator Gas specific Gravity (relative to air). If not defined, will use sg_g instead
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not defined, will use sg_sp instead
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point pressure (psia | barsa). If not provided, will attempt to calculate
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb | sm3/sm3)
        rsi: Initial solution GOR (scf/stb | sm3/sm3). If 0 (default), uses rsb
        zmethod: A string or z_method Enum class that specifies calculation method for gas Z-factor
        rsmethod: A string or rs_method Enum class that specifies calculation method for GOR
        cmethod: A string or c_method Enum class that specifies calculation method for gas critical properties
        denomethod: A string or deno_method Enum class that specifies calculation method for live oil density
        bomethod: A string or bo_method Enum class that specifies calculation method for oil FVF
        pbmethod: A string or pb_method Enum class that specifies calculation method for bubble point pressure
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3, rm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        if rsi > 0:
            rsi = rsi * SM3_PER_SM3_TO_SCF_PER_STB

    validate_pe_inputs(p=p, degf=degf)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

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

    if pb <= 0:
        pb = oil_pbub(
            api=api, degf=degf, rsb=rsb, sg_sp=sg_sp, pbmethod=pbmethod
        )
    if rsb <= 0:
        rsb = oil_rs_bub(
            api=api,
            degf=degf,
            pb=pb,
            sg_sp=sg_sp,
            rsmethod=rsmethod,
        )

    if rsi <= 0:
        rsi = rsb

    sg_o = oil_sg(api)

    rs = oil_rs(
        api=api, degf=degf, sg_sp=sg_sp, p=p,
        pb=pb, rsb=rsb, rsmethod=rsmethod, pbmethod=pbmethod,
    )

    bo = oil_bo(
        p=p, pb=pb, degf=degf, rs=rs, rsb=rsb,
        sg_sp=sg_sp, sg_g=sg_g, sg_o=sg_o,
        bomethod=bomethod, denomethod=denomethod,
    )

    if p >= pb:
        return bo  # No free gas above Pb

    # Bg in rcf/scf -> rb/scf
    bg = gas.gas_bg(p, sg_sp, degf, zmethod=zmethod, cmethod=cmethod) / CUFTperBBL

    bt = bo + (rsi - rs) * bg
    return bt
