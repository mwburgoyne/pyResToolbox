"""
Garcia (2001) density mixing rule for gas-dissolved aqueous solutions.

Single gas (Garcia Eq. 18):
    ρ = (1 + x2·M2/(M1·x1)) / (x2·V_phi/(M1·x1) + 1/ρ1)

Mixed gas (mole-fraction-weighted):
    V_phi_eff = Σ yi·V_phi_i
    M2_eff = Σ yi·M2_i
    where yi = mole fraction of gas i among dissolved gases (Σyi = 1)
    Then apply single-gas formula with V_phi_eff, M2_eff.

Brine support:
    For saline solutions, ρ1 is replaced with ρ_brine from Batzle & Wang (1992).
    V_phi values remain from the Plyasunov model (Garcia found salinity effects
    on V_phi are weak and within experimental uncertainty).

Units:
    T in K, P in MPa, S in weight fraction NaCl
    x2 = total dissolved gas mole fraction in liquid phase
    ρ in kg/m³, V_phi in cm³/mol, M in g/mol
"""

import numpy as np
from water_properties import rho_w, MW_WATER
from plyasunov_model import V_phi, gas_mw
from brine_properties import rho_brine


def density_single_gas(gas, x2, T, P, rho1=None, S=0.0):
    """
    Solution density with a single dissolved gas using Garcia Eq. 18.

    Parameters:
        gas: gas name string ('CO2', 'CH4', etc.)
        x2: mole fraction of dissolved gas in liquid phase
        T: temperature in K
        P: pressure in MPa
        rho1: solvent density in kg/m³ (computed from S if None)
        S: salinity as weight fraction NaCl (default 0 = pure water)

    Returns:
        solution density in kg/m³
    """
    if rho1 is None:
        rho1 = rho_brine(T, P, S) if S > 0 else rho_w(T, P)

    if x2 <= 0:
        return float(rho1)

    x1 = 1.0 - x2
    M1 = MW_WATER        # g/mol
    M2 = gas_mw(gas)     # g/mol
    vphi = V_phi(gas, T, P)  # cm³/mol

    # Garcia Eq. 18: ρ = (1 + x2·M2/(M1·x1)) / (x2·V_phi/(M1·x1) + 1/ρ1)
    # V_phi is in cm³/mol, M1 in g/mol, ρ1 in kg/m³
    # Need consistent units: convert V_phi to m³/mol? or ρ1 to g/cm³?
    #
    # Work in g/cm³ for clarity:
    # ρ1_gcc = ρ1 / 1000 [g/cm³]
    # Then: ρ [g/cm³] = (1 + x2·M2/(M1·x1)) / (x2·V_phi/(M1·x1) + 1/ρ1_gcc)
    # Convert result back to kg/m³ by multiplying by 1000.

    rho1_gcc = rho1 / 1000.0  # kg/m³ → g/cm³

    numerator = 1.0 + x2 * M2 / (M1 * x1)
    denominator = x2 * vphi / (M1 * x1) + 1.0 / rho1_gcc

    rho_gcc = numerator / denominator  # g/cm³
    return rho_gcc * 1000.0  # kg/m³


def density_mixed_gas(gas_dict, T, P, rho1=None, S=0.0):
    """
    Solution density with multiple dissolved gases.

    Uses mole-fraction-weighted effective V_phi and MW, then applies
    Garcia Eq. 18 with the effective values.

    Parameters:
        gas_dict: dict of {gas_name: x2_i} where x2_i is the mole fraction
                  of that gas in the liquid phase. Sum of all x2_i = total x2.
        T: temperature in K
        P: pressure in MPa
        rho1: solvent density in kg/m³ (computed from S if None)
        S: salinity as weight fraction NaCl (default 0 = pure water)

    Returns:
        solution density in kg/m³
    """
    if rho1 is None:
        rho1 = rho_brine(T, P, S) if S > 0 else rho_w(T, P)

    x2_total = sum(gas_dict.values())
    if x2_total <= 0:
        return float(rho1)

    # Compute mole-fraction-weighted effective properties
    # yi = x2_i / x2_total (fraction of gas i among all dissolved gases)
    vphi_eff = 0.0
    mw_eff = 0.0
    for gas, x2_i in gas_dict.items():
        yi = x2_i / x2_total
        vphi_eff += yi * V_phi(gas, T, P)
        mw_eff += yi * gas_mw(gas)

    # Apply Garcia Eq. 18 with effective values
    x1 = 1.0 - x2_total
    M1 = MW_WATER
    rho1_gcc = rho1 / 1000.0

    numerator = 1.0 + x2_total * mw_eff / (M1 * x1)
    denominator = x2_total * vphi_eff / (M1 * x1) + 1.0 / rho1_gcc

    rho_gcc = numerator / denominator
    return rho_gcc * 1000.0


# ============================================================================
# VISCOSITY CORRECTION
# ============================================================================
#
# Experimentally-calibrated viscosity corrections for dissolved gases.
# Each gas uses the best available experimental source:
#
#   CO2: Islam & Carlson (2012), Energy Fuels 26(8), 5330-5336
#        mu = mu_brine * (1 + 4.65 * x_CO2^1.0134)
#
#   CH4: Ostermann, Bloori & Dehghani (1985), SPE 14211
#        Confirmed by Ostermann et al. (1986), SPE 15081
#        mu_sat/mu_free = 1.109 - 5.98e-4*T + 1.0933e-6*T^2  (T in degF)
#        NOTE: SPE 14211 prints 1.0933e-5 (TYPO); verified 1.0933e-6 from
#        stated plateau values (1.060 at 100F, 1.044 at 150F, 1.028 at 250F).
#        Plateau above ~2000 psi; T-dependent, not mole-fraction-dependent.
#
#   H2S: Murphy & Gaines (1974), J. Chem. Eng. Data 19(4), 359-362
#        mu = mu_brine * (1 + 1.5 * x_H2S^1.0134)
#        Calibrated from 3-6% increase near saturation at 28-30 degC.
#        x_H2S at saturation ~ 0.027-0.033 (from S&W VLE at Murphy & Gaines
#        conditions). Implied 'a' ranges from 1.0 to 2.1 across the 5 data
#        points (large scatter; authors note the effect is "only slightly more
#        than experimental error"). Central estimate a=1.5.
#        Like CH4, effect is strongly T-dependent (nearly zero at 35 degC),
#        but only 5 data points spanning 28-35 degC — insufficient for a
#        T-dependent correlation. The constant a=1.5 is conservative at
#        higher T and approximate at lower T.
#
#   C2H6: NO EFFECT — Ostermann et al. (1986) SPE 15081
#   N2:   NO EFFECT — Murphy & Gaines (1974) control experiment
#   H2, C3H8, nC4H10: No data — no correction applied (conservative)
#
# NOTE: A previous density-based scaling approach (a_i = a_CO2 * drho%_i/drho%_CO2)
# was found to be fundamentally wrong. It predicted large viscosity DECREASES for
# CH4 (a=-11.0) and H2S (a=-0.9), but experiments show CH4 INCREASES viscosity
# by up to 6% and H2S INCREASES by 3-6%. The physical mechanism is hydrophobic
# hydration — dissolved gas molecules organize water into cage-like structures,
# increasing structural order and viscosity. This is unrelated to density effects.

_IC_A_CO2 = 4.65       # Islam-Carlson coefficient for CO2
_IC_A_H2S = 1.50       # Calibrated from Murphy & Gaines (1974), central estimate
_IC_B = 1.0134         # Islam-Carlson exponent (nearly linear)

# Ostermann (1985) Eq. 10 / (1986) Eq. 7 coefficients for CH4
# mu_sat/mu_free = c0 + c1*T + c2*T^2  where T in degF
# NOTE: SPE 14211 prints the quadratic coefficient as 1.0933e-5, but this
# is a TYPOGRAPHICAL ERROR. Back-calculation from the stated plateau values
# (1.060 at 100F, 1.044 at 150F, 1.028 at 250F) proves the correct
# exponent is 10^-6. SPE 15081 prints the linear term as -5.93e-4 (vs
# -5.98e-4 in SPE 14211); this minor difference does not affect results
# at the stated precision.
_OST_C0 = 1.109
_OST_C1 = -5.98e-4
_OST_C2 = 1.0933e-6


def _ostermann_ch4_plateau(degf):
    """
    CH4-saturated water viscosity ratio at plateau (P > ~2000 psi).

    Ostermann et al. (1985) SPE 14211, Eq. 10:
        mu_sat/mu_free = 1.109 - 5.98e-4*T + 1.0933e-5*T^2

    Validated values: 1.060 at 100 degF, 1.044 at 150 degF, 1.028 at 250 degF.
    Effect diminishes with temperature but remains positive (>1.0) across
    the full reservoir temperature range.

    Parameters:
        degf: temperature in degrees Fahrenheit

    Returns:
        viscosity ratio mu_sat/mu_free (always >= 1.0)
    """
    ratio = _OST_C0 + _OST_C1 * degf + _OST_C2 * degf ** 2
    return max(ratio, 1.0)


def viscosity_correction_single(gas, x2, degf=None):
    """
    Viscosity multiplier for a single dissolved gas.

    Uses experimentally-calibrated corrections:
      - CO2: Islam-Carlson (2012) power-law
      - CH4: Ostermann (1985) T-dependent plateau (requires degf)
      - H2S: Calibrated from Murphy & Gaines (1974)
      - All others: no correction (1.0)

    mu_corrected = mu_brine * viscosity_correction_single(gas, x2, degf)

    Parameters:
        gas: gas name string ('CO2', 'CH4', 'H2S', etc.)
        x2: mole fraction of dissolved gas in liquid phase
        degf: temperature in degrees Fahrenheit (required for CH4)

    Returns:
        viscosity multiplier (>= 1.0 for all calibrated gases)
    """
    if x2 <= 0:
        return 1.0

    gas = gas.upper()

    if gas == 'CO2':
        return 1.0 + _IC_A_CO2 * x2 ** _IC_B

    if gas == 'H2S':
        return 1.0 + _IC_A_H2S * x2 ** _IC_B

    if gas == 'CH4':
        if degf is None:
            raise ValueError("degf (temperature) is required for CH4 viscosity correction")
        # Ostermann plateau value — the full effect at saturation.
        # The plateau is reached above ~2000 psi CH4 partial pressure.
        # For typical reservoir conditions (P > 2000 psi), this is appropriate.
        # For undersaturated conditions, the effect could be scaled by
        # Rsw/Rsw_sat, but Ostermann's data shows the plateau is reached
        # at relatively low pressures, so the full correction is used here.
        return _ostermann_ch4_plateau(degf)

    # C2H6, N2: experimentally confirmed no effect (Ostermann 1986, Murphy & Gaines 1974)
    # H2, C3H8, NC4H10: no data available — conservative (no correction)
    return 1.0


def viscosity_correction_mixed(gas_dict, degf=None):
    """
    Viscosity multiplier for multiple dissolved gases.

    Applies multiplicative corrections:
        mu = mu_brine * Product_i(correction_i)

    Parameters:
        gas_dict: dict of {gas_name: x2_i}
        degf: temperature in degrees Fahrenheit (required if CH4 present)

    Returns:
        combined viscosity multiplier
    """
    factor = 1.0
    for gas, x2 in gas_dict.items():
        if x2 > 0:
            factor *= viscosity_correction_single(gas, x2, degf=degf)
    return factor


def density_change_pct(gas_or_dict, x2_or_none, T, P, S=0.0):
    """
    Percentage density change relative to solvent (water or brine).

    Parameters:
        gas_or_dict: either a gas name string, or a dict {gas: x2_i}
        x2_or_none: mole fraction if gas_or_dict is a string, None if dict
        T, P: temperature (K) and pressure (MPa)
        S: salinity as weight fraction NaCl (default 0 = pure water)

    Returns:
        (rho_solution, rho_solvent, percent_change)
    """
    rho1 = rho_brine(T, P, S) if S > 0 else rho_w(T, P)
    if isinstance(gas_or_dict, dict):
        rho_sol = density_mixed_gas(gas_dict=gas_or_dict, T=T, P=P, rho1=rho1)
    else:
        rho_sol = density_single_gas(gas_or_dict, x2_or_none, T, P, rho1=rho1)
    pct = 100.0 * (rho_sol - rho1) / rho1
    return rho_sol, rho1, pct


# ============================================================================
# MAIN - VALIDATION
# ============================================================================

if __name__ == "__main__":
    # Validate against Garcia (2001) CO2 results
    # Garcia reports ~2.5% density increase at x_CO2 = 0.05
    print("=== Single Gas Density Tests ===\n")

    T_test = 298.15  # 25°C
    P_test = 10.0    # 10 MPa (representative reservoir P)

    rho1 = rho_w(T_test, P_test)
    print(f"Pure water at T={T_test}K, P={P_test}MPa: ρ = {rho1:.4f} kg/m³")
    print()

    gases = ['CO2', 'CH4', 'H2S', 'N2', 'H2', 'C2H6', 'C3H8', 'NC4H10']
    x2_test = 0.02  # 2 mol% dissolved gas

    print(f"Dissolved gas mole fraction x2 = {x2_test}")
    print(f"{'Gas':<8} {'MW':>6} {'V_phi':>8} {'ρ_sol':>10} {'Δρ':>8} {'%change':>8}")
    print("-" * 55)

    for gas in gases:
        vphi = V_phi(gas, T_test, P_test)
        rho_sol = density_single_gas(gas, x2_test, T_test, P_test, rho1=rho1)
        delta = rho_sol - rho1
        pct = 100 * delta / rho1
        mw = gas_mw(gas)
        print(f"{gas:<8} {mw:6.2f} {vphi:8.2f} {rho_sol:10.4f} {delta:+8.4f} {pct:+7.3f}%")

    # Test: CO2 at x2=0.05 should give ~2.5% increase (Garcia Fig. 3)
    print(f"\n--- CO2 at x2=0.05 (Garcia reports ~2.5% increase) ---")
    rho_co2, _, pct_co2 = density_change_pct('CO2', 0.05, T_test, P_test)
    print(f"  ρ = {rho_co2:.4f} kg/m³, Δρ/ρ = {pct_co2:+.3f}%")

    # Mixed gas test
    print("\n=== Mixed Gas Test ===\n")
    gas_mix = {'CO2': 0.01, 'CH4': 0.005, 'H2S': 0.003}
    rho_mix = density_mixed_gas(gas_mix, T_test, P_test, rho1=rho1)
    delta_mix = rho_mix - rho1
    pct_mix = 100 * delta_mix / rho1
    total_x2 = sum(gas_mix.values())
    print(f"Gas mix: {gas_mix}")
    print(f"Total x2 = {total_x2}")
    print(f"ρ_mix = {rho_mix:.4f} kg/m³, Δρ = {delta_mix:+.4f}, %change = {pct_mix:+.3f}%")

    # Brine tests
    print("\n=== Brine + Gas Tests ===\n")
    S_test = 0.10  # 10 wt% NaCl
    rho_brine_val = rho_brine(T_test, P_test, S_test)
    print(f"Brine at T={T_test}K, P={P_test}MPa, S={S_test*100}wt%: ρ = {rho_brine_val:.4f} kg/m³")
    print(f"\nDensity effect of dissolved gas in brine vs pure water (x2=0.02):")
    print(f"{'Gas':<8} {'ρ(water)':>10} {'Δρ%(w)':>8} {'ρ(brine)':>10} {'Δρ%(b)':>8}")
    print("-" * 50)
    for gas in ['CO2', 'CH4', 'H2S', 'N2']:
        rho_w_sol, rho_w_base, pct_w = density_change_pct(gas, x2_test, T_test, P_test, S=0.0)
        rho_b_sol, rho_b_base, pct_b = density_change_pct(gas, x2_test, T_test, P_test, S=S_test)
        print(f"{gas:<8} {rho_w_sol:10.4f} {pct_w:+7.3f}% {rho_b_sol:10.4f} {pct_b:+7.3f}%")

    # Temperature sweep for CO2
    print("\n=== CO2 Density vs Temperature (x2=0.02, P=30 MPa) ===\n")
    print(f"{'T(°C)':<8} {'ρ_water':>10} {'ρ_CO2aq':>10} {'Δρ%':>8}")
    print("-" * 40)
    for T_C in [25, 50, 100, 150, 200, 250]:
        T_K = T_C + 273.15
        rho_w_val = rho_w(T_K, 30.0)
        rho_sol_val = density_single_gas('CO2', 0.02, T_K, 30.0, rho1=rho_w_val)
        pct_val = 100 * (rho_sol_val - rho_w_val) / rho_w_val
        print(f"{T_C:<8} {rho_w_val:10.4f} {rho_sol_val:10.4f} {pct_val:+7.3f}%")

    # ================================================================
    # VISCOSITY VALIDATION
    # ================================================================
    print("\n=== Viscosity Correction Tests ===\n")

    # CO2: Islam-Carlson (2012) at x=0.02
    co2_vis = viscosity_correction_single('CO2', 0.02)
    print(f"CO2 at x=0.02: factor = {co2_vis:.6f} (expected ~1.09)")

    # H2S: Murphy & Gaines calibrated (a=0.8)
    h2s_vis_005 = viscosity_correction_single('H2S', 0.05)
    h2s_vis_020 = viscosity_correction_single('H2S', 0.20)
    print(f"H2S at x=0.05: factor = {h2s_vis_005:.6f} (expect ~+4% increase)")
    print(f"H2S at x=0.20: factor = {h2s_vis_020:.6f} (expect ~+15.5%)")

    # CH4: Ostermann plateau values
    print("\nCH4 Ostermann plateau vs temperature:")
    print(f"  {'T(degF)':<10} {'Plateau':>10} {'Expected':>10}")
    print(f"  {'-'*32}")
    for degf, expected in [(100, 1.060), (150, 1.044), (250, 1.028)]:
        plateau = _ostermann_ch4_plateau(degf)
        print(f"  {degf:<10} {plateau:10.4f} {expected:10.3f}")

    # Gases with no effect
    print(f"\nC2H6 at x=0.02: factor = {viscosity_correction_single('C2H6', 0.02):.4f} (expect 1.0)")
    print(f"N2 at x=0.02:   factor = {viscosity_correction_single('N2', 0.02):.4f} (expect 1.0)")
    print(f"H2 at x=0.02:   factor = {viscosity_correction_single('H2', 0.02):.4f} (expect 1.0)")

    # Mixed gas test
    gas_mix_v = {'CO2': 0.01, 'CH4': 0.005, 'H2S': 0.003}
    mixed_vis = viscosity_correction_mixed(gas_mix_v, degf=150)
    print(f"\nMixed viscosity ({gas_mix_v}, T=150F): factor = {mixed_vis:.6f}")
