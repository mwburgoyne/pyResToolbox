"""Gas hydrate formation prediction, water balance, and inhibitor calculations.

Carved out of ``gas.py`` (v3.5.0). Public API is unchanged — ``gas_hydrate`` and
``HydrateResult`` are re-exported from ``pyrestoolbox.gas`` via ``gas.py``.

The only intra-``gas`` dependency is ``gas_water_content``; it is imported lazily
inside ``gas_hydrate`` to avoid a circular import with ``gas.gas``.
"""

import math
from dataclasses import dataclass

from pyrestoolbox.classes import hyd_method, inhibitor
from pyrestoolbox.shared_fns import validate_pe_inputs
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.constants import (
    BAR_TO_PSI, PSI_TO_BAR, degc_to_degf, degf_to_degc,
    STB_PER_MMSCF_TO_SM3_PER_SM3, SM3_PER_SM3_TO_STB_PER_MMSCF,
    LB_PER_MMSCF_TO_KG_PER_SM3, GAL_PER_MMSCF_TO_L_PER_SM3)

# --- Motiee (1991) hydrate formation temperature ---
# Hydrocarbon Processing 70, pp 98-99. degF-output, psia-input form:
# T(degF) = b0 + b1*log10(P) + b2*log10(P)^2 + b3*g + b4*g^2 + b5*g*log10(P)
# Coefficient set per the degF/psia statement of the correlation (intercept
# -238.24469); cross-checked against GPSA chart values (1000 psia, sg 0.6
# gives 56.3 degF vs chart ~60-62 degF, consistent with Motiee's documented
# 3-7 degF under-prediction of the GPSA chart). Other unit/intercept variants
# circulating in the literature (degC/kPa or degF/psia with -283.24469) fail
# this sanity check by 30-50 degF.
_MOTIEE = (-238.24469, 78.99181, -5.352544, 349.47324, -150.85396, -27.604065)

# --- Towler & Mokhatab (2005) hydrate formation temperature ---
# Hydrocarbon Processing 84, pp 61-62
_TOWLER = (13.47, 34.27, -1.675, -20.35)

# --- Hydrate formation pressure search bounds ---
_HFP_P_LO = 14.696   # ~1 atm (psia)
_HFP_P_HI = 15000.0   # Upper search bound (psia)

# Østergaard et al. (2005) coefficients: delta_T(degC) = C1*w + C2*w^2 + C3*w^3, w in wt% (0-100)
_OSTERGAARD = {
    inhibitor.MEOH: (0.4411, -0.0033, 6.476e-5),
    inhibitor.MEG:  (0.2533, -0.0009, 2.222e-5),
    inhibitor.DEG:  (0.1833, -0.0004, 6.667e-6),
    inhibitor.TEG:  (0.1333, -0.0002, 3.333e-6),
    inhibitor.ETOH: (0.3750, -0.0020, 3.500e-5),
}


# Inhibitor physical density at 20 degC (g/cm3) — for volumetric injection rate
_INHIBITOR_DENSITY = {
    inhibitor.MEOH: 0.791,
    inhibitor.MEG:  1.110,
    inhibitor.DEG:  1.117,
    inhibitor.TEG:  1.125,
    inhibitor.ETOH: 0.789,
}

# Conversion: g/cm3 to lb/gal
_GCM3_TO_LB_PER_GAL = 8.34540445

# Water mass at standard conditions: lb per stb
_WATER_LB_PER_STB = 350.2

# Maximum valid wt% for each inhibitor (aqueous phase concentration)
_MAX_WT_PCT = {
    inhibitor.MEOH: 25.0,
    inhibitor.MEG:  70.0,
    inhibitor.DEG:  70.0,
    inhibitor.TEG:  50.0,
    inhibitor.ETOH: 30.0,
}


@dataclass
class HydrateResult:
    """Result of gas hydrate formation prediction.

    Hydrate assessment (HFT, HFP, subcooling, inhibitor) is evaluated at the
    operating point (p, degf). Water balance is evaluated between reservoir
    conditions (p_res, degf_res) and the operating point.

    Attributes — Hydrate Assessment
    --------------------------------
    hft : float
        Hydrate formation temperature at operating pressure (degF | degC).
    hfp : float
        Hydrate formation pressure at operating temperature (psia | barsa).
    subcooling : float
        HFT - T_operating (degF | degC delta). Positive means inside hydrate window.
    in_hydrate_window : bool
        True if operating temperature is below the hydrate formation temperature.
    inhibited_hft : float
        HFT after inhibitor depression (degF | degC), or NaN if no inhibitor.
    inhibitor_depression : float
        Temperature depression from inhibitor (degF | degC delta), or 0.
    required_inhibitor_wt_pct : float
        Wt% inhibitor in aqueous phase (water + inhibitor) needed to bring HFT
        below operating temperature, or 0. Capped at the physical maximum for
        the selected inhibitor type.
    max_inhibitor_wt_pct : float
        Maximum valid wt% for the selected inhibitor type (MEOH: 25%, MEG: 70%,
        DEG: 70%, TEG: 50%, ETOH: 30%). 0 if no inhibitor specified.
    inhibitor_underdosed : bool
        True if required_inhibitor_wt_pct exceeds max_inhibitor_wt_pct,
        meaning this inhibitor type cannot provide sufficient depression
        even at its physical maximum concentration. Does NOT indicate
        whether the applied inhibitor_wt_pct is sufficient — compare
        inhibited_hft to operating temperature to check applied-dose
        protection.

    Attributes — Water Balance
    --------------------------
    The gas leaves the reservoir saturated with vaporized water at reservoir
    P,T. At the operating point (lower P,T), the gas can hold less water vapor,
    so the excess condenses as liquid. Free water is any additional liquid water
    entrained in the gas stream from the reservoir (user-specified).

    water_vaporized_res : float
        Equilibrium vaporized water in the gas at reservoir P,T
        (stb/MMscf | sm3/sm3). When p_res/degf_res not provided, equals
        water_vaporized_op (no temperature/pressure change, no condensation).
    water_vaporized_op : float
        Equilibrium vaporized water in the gas at operating P,T
        (stb/MMscf | sm3/sm3).
    water_condensed : float
        Water condensed from vapor between reservoir and operating conditions
        = max(water_vaporized_res - water_vaporized_op, 0) (stb/MMscf | sm3/sm3).
    free_water : float
        Free liquid water influx from reservoir (= additional_water input)
        (stb/MMscf | sm3/sm3).
    total_liquid_water : float
        Total liquid water at operating point = water_condensed + free_water
        (stb/MMscf | sm3/sm3). This is the water that must be treated with
        inhibitor to prevent hydrate formation.

    Attributes — Inhibitor Injection Rate
    --------------------------------------
    Injection rate is based on the total liquid water at the operating point
    and the required inhibitor concentration (wt% in the aqueous phase).

    inhibitor_mass_rate : float
        Required inhibitor mass injection rate (lb/MMscf | kg/sm3 if metric).
        Zero when outside hydrate window, no inhibitor, or no liquid water.
    inhibitor_vol_rate : float
        Required inhibitor volume injection rate (gal/MMscf | L/sm3 if metric).
        Zero when outside hydrate window, no inhibitor, or no liquid water.
    """
    hft: float
    hfp: float
    subcooling: float
    in_hydrate_window: bool
    inhibited_hft: float
    inhibitor_depression: float
    required_inhibitor_wt_pct: float
    max_inhibitor_wt_pct: float
    inhibitor_underdosed: bool
    water_vaporized_res: float
    water_vaporized_op: float
    water_condensed: float
    free_water: float
    total_liquid_water: float
    inhibitor_mass_rate: float
    inhibitor_vol_rate: float


def _motiee_hft(p_psia, sg):
    """Motiee (1991) hydrate formation temperature.

    Published form gives T directly in degF from log10 of pressure in psia:
    T(degF) = -238.24469 + 78.99181*log10(P_psia) - 5.352544*log10(P_psia)^2
              + 349.47324*gamma - 150.85396*gamma^2 - 27.604065*gamma*log10(P_psia)

    Valid for 100-4000 psia and gas SG 0.55-0.90. Fitted to the GPSA (Katz)
    gravity chart; tends to under-predict the chart by roughly 3-7 degF.

    Reference: Motiee, M. (1991). Hydrocarbon Processing 70, pp 98-99.
    """
    log_p = math.log10(p_psia)
    return (_MOTIEE[0]
            + _MOTIEE[1] * log_p
            + _MOTIEE[2] * log_p * log_p
            + _MOTIEE[3] * sg
            + _MOTIEE[4] * sg * sg
            + _MOTIEE[5] * sg * log_p)


def _towler_mokhatab_hft(p_psia, sg):
    """Towler & Mokhatab (2005) hydrate formation temperature.

    T(degF) = 13.47*ln(P_psia) + 34.27*ln(gamma) - 1.675*ln(P_psia)*ln(gamma) - 20.35

    Reference: Towler, B.F. & Mokhatab, S. (2005). Hydrocarbon Processing 84, pp 61-62.
    """
    ln_p = math.log(p_psia)
    ln_sg = math.log(sg)
    return _TOWLER[0] * ln_p + _TOWLER[1] * ln_sg + _TOWLER[2] * ln_p * ln_sg + _TOWLER[3]


def _hydrate_formation_press(degf_target, sg, hft_fn):
    """Invert the supplied HFT correlation (hft_fn) to find the hydrate
    formation pressure via bisection.

    Returns pressure in psia, or NaN if target T is outside correlation range.
    """
    p_lo = _HFP_P_LO
    p_hi = _HFP_P_HI

    t_lo = hft_fn(p_lo, sg)
    t_hi = hft_fn(p_hi, sg)

    # Check that solution exists within bounds
    if degf_target < t_lo or degf_target > t_hi:
        return float('nan')

    for _ in range(100):
        p_mid = (p_lo + p_hi) / 2.0
        t_mid = hft_fn(p_mid, sg)
        if abs(t_mid - degf_target) < 0.001:
            break
        if t_mid < degf_target:
            p_lo = p_mid
        else:
            p_hi = p_mid

    return (p_lo + p_hi) / 2.0


def _ostergaard_depression(wt_pct, inh):
    """Østergaard et al. (2005) temperature depression in degC.

    delta_T(degC) = C1*w + C2*w^2 + C3*w^3, where w = wt% (0-100 scale).

    Reference: Østergaard, K.K. et al. (2005). J. Pet. Sci. Eng. 48, pp 70-80.
    """
    c1, c2, c3 = _OSTERGAARD[inh]
    return c1 * wt_pct + c2 * wt_pct**2 + c3 * wt_pct**3


def _required_concentration(depression_degc, inh):
    """Newton-Raphson inversion of Østergaard cubic to find required wt%.

    Returns wt% (0-100 scale). Returns 0 if depression <= 0.
    """
    if depression_degc <= 0:
        return 0.0

    c1, c2, c3 = _OSTERGAARD[inh]

    # Initial guess from linear term
    w = depression_degc / c1 if c1 > 0 else 20.0

    for _ in range(50):
        f = c1 * w + c2 * w**2 + c3 * w**3 - depression_degc
        fp = c1 + 2 * c2 * w + 3 * c3 * w**2
        if abs(fp) < 1e-15:
            break
        w_new = w - f / fp
        if w_new < 0:
            w_new = w / 2.0
        if abs(w_new - w) < 1e-6:
            w = w_new
            break
        w = w_new

    return max(float(w), 0.0)


def gas_hydrate(
    p: float,
    degf: float,
    sg: float,
    hydmethod: str = 'TOWLER',
    inhibitor_type: str = None,
    inhibitor_wt_pct: float = 0,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    p_res: float = None,
    degf_res: float = None,
    additional_water: float = 0,
    metric: bool = False,
) -> HydrateResult:
    """ Returns gas hydrate formation prediction, water balance, and inhibitor calculations.

        Hydrate assessment (HFT, HFP, subcooling, inhibitor) is evaluated at
        the operating point (p, degf) — typically the wellhead or coldest point
        in the production system. Water balance is evaluated between reservoir
        conditions (p_res, degf_res) and the operating point to determine how
        much water condenses from vapor, how much was always liquid (free water),
        and the total liquid water that must be treated with inhibitor.

        p: Operating pressure at hydrate assessment point, e.g. wellhead
           (psia | barsa if metric=True)
        degf: Operating temperature at hydrate assessment point
              (degF | degC if metric=True)
        sg: Gas specific gravity (air = 1.0)
        hydmethod: Hydrate formation correlation.
                   'TOWLER': Towler & Mokhatab (2005)
                   'MOTIEE': Motiee (1991)
                   Defaults to 'TOWLER'
        inhibitor_type: Thermodynamic hydrate inhibitor type (optional).
                        'MEOH' (Methanol), 'MEG' (Monoethylene Glycol),
                        'DEG' (Diethylene Glycol), 'TEG' (Triethylene Glycol),
                        'ETOH' (Ethanol). None = no inhibitor
        inhibitor_wt_pct: Weight percent of inhibitor in aqueous phase (water
                          + inhibitor, 0-100). Defaults to 0
        co2: CO2 mole fraction (0-1). For composition-aware water content via
             SoreideWhitson. Defaults to 0
        h2s: H2S mole fraction (0-1). Defaults to 0
        n2: N2 mole fraction (0-1). Defaults to 0
        h2: H2 mole fraction (0-1). Defaults to 0
        p_res: Reservoir pressure where gas was last in equilibrium with water
               (psia | barsa if metric=True). Determines how much water the gas
               carries as vapor from the reservoir. If None, uses p (operating
               pressure — no condensation). Defaults to None
        degf_res: Reservoir temperature where gas was last in equilibrium with
                  water (degF | degC if metric=True). Determines how much water
                  the gas carries as vapor from the reservoir. If None, uses degf
                  (operating temperature — no condensation). Defaults to None
        additional_water: Free liquid water entrained in the gas stream from the
                          reservoir, e.g. from mobile formation water
                          (stb/MMscf | sm3/sm3 if metric). This water was never
                          vaporized — it travels with the gas as liquid. Added to
                          condensed water for inhibitor dosing. Defaults to 0
        metric: If True, input/output in Eclipse METRIC units (barsa, degC).
                Defaults to False (FIELD: psia, degF)

    Returns a HydrateResult dataclass. See HydrateResult docstring for full
    field descriptions.
    """
    # Composition-aware water content uses gas_water_content; imported lazily to
    # avoid a circular import (gas.gas imports this module at top level).
    from pyrestoolbox.gas.gas import gas_water_content

    # Convert metric inputs to oilfield
    if metric:
        p_psia = p * BAR_TO_PSI
        degf_of = degc_to_degf(degf)
        additional_water_stb = additional_water * SM3_PER_SM3_TO_STB_PER_MMSCF
        p_res_psia = p_res * BAR_TO_PSI if p_res is not None else None
        degf_res_of = degc_to_degf(degf_res) if degf_res is not None else None
    else:
        p_psia = p
        degf_of = degf
        additional_water_stb = additional_water
        p_res_psia = p_res
        degf_res_of = degf_res

    # Water content P,T: use reservoir conditions if provided, else operating
    p_wc = p_res_psia if p_res_psia is not None else p_psia
    degf_wc = degf_res_of if degf_res_of is not None else degf_of

    # Validate inputs
    validate_pe_inputs(p=p_psia, degf=degf_of, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    if p_wc != p_psia:
        validate_pe_inputs(p=p_wc)
    if inhibitor_wt_pct < 0 or inhibitor_wt_pct >= 100:
        raise ValueError("Inhibitor wt% must be in range [0, 100)")
    if additional_water < 0:
        raise ValueError("additional_water must be >= 0")

    # Resolve hydrate method
    hydmethod = validate_methods(["hydmethod"], [hydmethod])

    # Select HFT function
    if hydmethod == hyd_method.MOTIEE:
        hft_fn = _motiee_hft
    else:
        hft_fn = _towler_mokhatab_hft

    # Compute HFT and HFP
    hft_degf = hft_fn(p_psia, sg)
    hfp_psia = _hydrate_formation_press(degf_of, sg, hft_fn)

    # Subcooling and hydrate window
    subcooling_degf = hft_degf - degf_of
    in_window = degf_of < hft_degf

    # --- Water balance ---
    # Compute vaporized water at both reservoir and operating conditions.
    # Condensed water = what dropped out of vapor between the two points.
    # Free water = liquid water entrained from reservoir (user input).
    # Total liquid = condensed + free = what needs inhibitor treatment.
    has_composition = co2 > 0 or h2s > 0 or n2 > 0 or h2 > 0

    def _water_content_at(p_eval, degf_eval):
        """Compute equilibrium vaporized water content at given P,T (stb/MMscf)."""
        if has_composition:
            from pyrestoolbox.brine import SoreideWhitson
            sw = SoreideWhitson(
                pres=p_eval, temp=degf_eval, ppm=0,
                y_CO2=co2, y_H2S=h2s, y_N2=n2, y_H2=h2,
                sg=sg, metric=False,
            )
            return float(sw.water_content['stb_mmscf'])
        else:
            return float(gas_water_content(p=p_eval, degf=degf_eval, salinity=0, metric=False))

    # Vaporized water at operating point (always needed)
    wc_op = _water_content_at(p_psia, degf_of)

    # Vaporized water at reservoir (= operating if no reservoir P,T given)
    if p_res_psia is not None or degf_res_of is not None:
        wc_res = _water_content_at(p_wc, degf_wc)
    else:
        wc_res = wc_op  # no reservoir specified → no condensation

    # Condensed water: what dropped out of the gas between reservoir and operating
    condensed = max(wc_res - wc_op, 0.0)

    # Free water: liquid water entrained from reservoir (user input)
    free_water_stb = additional_water_stb

    # Total liquid water at operating point: this is what needs inhibitor
    total_liquid = condensed + free_water_stb

    # --- Inhibitor calculations ---
    inhibited_hft_degf = float('nan')
    depression_degf = 0.0
    required_wt_pct = 0.0
    max_wt_pct = 0.0
    underdosed = False
    inh_mass_rate = 0.0
    inh_vol_rate = 0.0

    if inhibitor_type is not None:
        # Resolve inhibitor enum
        inh = validate_methods(["inhibitor"], [inhibitor_type])
        max_wt_pct = _MAX_WT_PCT[inh]

        if inhibitor_wt_pct > 0:
            # Østergaard depression (in degC), convert to degF delta
            depression_degc = _ostergaard_depression(inhibitor_wt_pct, inh)
            depression_degf = depression_degc * 9.0 / 5.0
            inhibited_hft_degf = hft_degf - depression_degf
        else:
            inhibited_hft_degf = hft_degf
            depression_degf = 0.0

        # Required concentration to bring HFT below operating T (with capping)
        if degf_of < hft_degf:
            needed_depression_degc = (hft_degf - degf_of) * 5.0 / 9.0
            raw_wt_pct = _required_concentration(needed_depression_degc, inh)
            if raw_wt_pct > max_wt_pct:
                required_wt_pct = max_wt_pct
                underdosed = True
            else:
                required_wt_pct = raw_wt_pct
                underdosed = False
        else:
            required_wt_pct = 0.0

        # --- Injection rate ---
        # Inhibitor treats the total liquid water at operating conditions
        if in_window and required_wt_pct > 0 and total_liquid > 0:
            w_frac = required_wt_pct / 100.0
            liquid_mass_lb = total_liquid * _WATER_LB_PER_STB  # lb per MMscf gas
            inh_mass_rate = liquid_mass_lb * w_frac / (1.0 - w_frac)  # lb/MMscf
            density_lb_per_gal = _INHIBITOR_DENSITY[inh] * _GCM3_TO_LB_PER_GAL
            inh_vol_rate = inh_mass_rate / density_lb_per_gal  # gal/MMscf

    # Convert outputs if metric
    _wc_conv = STB_PER_MMSCF_TO_SM3_PER_SM3
    if metric:
        hft_out = degf_to_degc(hft_degf)
        hfp_out = hfp_psia * PSI_TO_BAR if not math.isnan(hfp_psia) else float('nan')
        subcooling_out = subcooling_degf * 5.0 / 9.0
        depression_out = depression_degf * 5.0 / 9.0
        inhibited_hft_out = degf_to_degc(inhibited_hft_degf) if not math.isnan(inhibited_hft_degf) else float('nan')
        wc_res_out = wc_res * _wc_conv
        wc_op_out = wc_op * _wc_conv
        condensed_out = condensed * _wc_conv
        free_water_out = free_water_stb * _wc_conv
        total_liquid_out = total_liquid * _wc_conv
        inh_mass_rate_out = inh_mass_rate * LB_PER_MMSCF_TO_KG_PER_SM3
        inh_vol_rate_out = inh_vol_rate * GAL_PER_MMSCF_TO_L_PER_SM3
    else:
        hft_out = hft_degf
        hfp_out = hfp_psia
        subcooling_out = subcooling_degf
        depression_out = depression_degf
        inhibited_hft_out = inhibited_hft_degf
        wc_res_out = wc_res
        wc_op_out = wc_op
        condensed_out = condensed
        free_water_out = free_water_stb
        total_liquid_out = total_liquid
        inh_mass_rate_out = inh_mass_rate
        inh_vol_rate_out = inh_vol_rate

    return HydrateResult(
        hft=hft_out,
        hfp=hfp_out,
        subcooling=subcooling_out,
        in_hydrate_window=in_window,
        inhibited_hft=inhibited_hft_out,
        inhibitor_depression=depression_out,
        required_inhibitor_wt_pct=required_wt_pct,
        max_inhibitor_wt_pct=max_wt_pct,
        inhibitor_underdosed=underdosed,
        water_vaporized_res=wc_res_out,
        water_vaporized_op=wc_op_out,
        water_condensed=condensed_out,
        free_water=free_water_out,
        total_liquid_water=total_liquid_out,
        inhibitor_mass_rate=inh_mass_rate_out,
        inhibitor_vol_rate=inh_vol_rate_out,
    )
