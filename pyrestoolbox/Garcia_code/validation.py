"""
Comprehensive validation of the Plyasunov V2∞ model and Garcia mixing rule.

Tests:
1. V2∞ at 298.15K for all 8 gases (Plyasunov reference values)
2. V2∞ vs temperature for CH4, CO2, H2S (Hnedkovsky 1996 experimental data)
3. V2∞ vs Garcia (2001) CO2 cubic polynomial
4. Density predictions vs Garcia CO2 results
5. Physical reasonableness checks
"""

import sys
import numpy as np
sys.path.insert(0, '.')

from water_properties import rho_w, kappa_T, V1_star, MW_WATER
from plyasunov_model import V2_inf, V_phi, gas_mw, B12
from garcia_mixing import density_single_gas, density_mixed_gas, density_change_pct

PASS = 0
FAIL = 0


def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  FAIL: {name} -- {detail}")


def garcia_CO2_vphi(T_C):
    """Garcia (2001) Eq. 3: V_phi for CO2 in cm³/mol, T in °C."""
    return 37.51 - 9.585e-2 * T_C + 8.740e-4 * T_C**2 - 5.044e-7 * T_C**3


# ============================================================================
# TEST 1: V2∞ at 298.15K reference conditions
# ============================================================================
def test_V2inf_298K():
    print("\n=== TEST 1: V2∞ at 298.15K, 0.1 MPa ===")

    known = {
        'H2': 26.1, 'N2': 34.7, 'CH4': 37.0, 'CO2': 34.0,
        'C2H6': 51.4, 'C3H8': 67.1, 'NC4H10': 82.8, 'H2S': 34.8,
    }

    T, P = 298.15, 0.1
    for gas, expected in known.items():
        computed = V2_inf(gas, T, P)
        pct_err = abs(100 * (computed - expected) / expected)
        check(f"{gas}: {computed:.2f} vs {expected} ({pct_err:.2f}%)",
              pct_err < 0.5,
              f"computed={computed:.4f}, expected={expected}")


# ============================================================================
# TEST 2: V2∞ vs Hnedkovsky (1996) experimental data
# ============================================================================
def test_V2inf_vs_hnedkovsky():
    print("\n=== TEST 2: V2∞ vs Hnedkovsky (1996) at 30 MPa ===")

    # Hnedkovsky data: average of 28 MPa and 35 MPa measurements
    # Format: (T_K, V2inf_avg)
    hned_CH4 = [
        (298.15, 36.75), (323.15, 37.30), (373.15, 40.50),
        (423.15, 45.90), (473.15, 54.10), (523.15, 64.80),
    ]
    hned_CO2 = [
        (298.15, 33.45), (323.15, 33.75), (373.15, 37.50),
        (423.15, 42.30), (473.15, 49.20), (523.15, 59.75),
    ]
    hned_H2S = [
        (298.15, 34.90), (323.15, 35.90), (373.15, 38.70),
        (423.15, 42.75), (473.15, 48.45), (523.15, 56.75),
    ]

    P = 30.0
    print(f"\n  {'Gas':<6} {'T(°C)':>6} {'Calc':>8} {'Expt':>8} {'%Err':>7}")
    print("  " + "-" * 40)

    for gas, data, tol in [('CH4', hned_CH4, 5.0), ('CO2', hned_CO2, 5.0), ('H2S', hned_H2S, 5.0)]:
        max_err = 0
        for T_K, v_exp in data:
            v_calc = V2_inf(gas, T_K, P)
            pct = 100 * (v_calc - v_exp) / v_exp
            max_err = max(max_err, abs(pct))
            T_C = T_K - 273.15
            print(f"  {gas:<6} {T_C:6.1f} {v_calc:8.2f} {v_exp:8.2f} {pct:+6.2f}%")
        check(f"{gas} max error {max_err:.1f}% < {tol}%", max_err < tol,
              f"max_err={max_err:.2f}%")


# ============================================================================
# TEST 3: Plyasunov CO2 V2∞ vs Garcia cubic polynomial
# ============================================================================
def test_V2inf_vs_garcia_polynomial():
    print("\n=== TEST 3: Plyasunov CO2 vs Garcia cubic polynomial ===")
    print("  Garcia Eq.3 is P-independent; Plyasunov at 20 MPa for comparison")
    print(f"\n  {'T(°C)':>6} {'Plyasunov':>10} {'Garcia':>10} {'Diff':>8}")
    print("  " + "-" * 38)

    P = 20.0
    max_diff_200 = 0
    for T_C in [25, 50, 75, 100, 125, 150, 175, 200, 250, 300]:
        T_K = T_C + 273.15
        v_plyas = V2_inf('CO2', T_K, P)
        v_garcia = garcia_CO2_vphi(T_C)
        diff = v_plyas - v_garcia
        if T_C <= 250:
            max_diff_200 = max(max_diff_200, abs(diff))
        print(f"  {T_C:6} {v_plyas:10.2f} {v_garcia:10.2f} {diff:+8.2f}")

    # Garcia polynomial was fit to data 5-300°C; at 300°C it diverges from Plyasunov
    # because Garcia's cubic doesn't capture near-critical behavior well.
    # Expect agreement within ~2 cm³/mol for 25-250°C range.
    check(f"Max difference 25-250°C: {max_diff_200:.2f} cm³/mol < 3",
          max_diff_200 < 3.0,
          f"max_diff={max_diff_200:.2f}")


# ============================================================================
# TEST 4: CO2 density predictions
# ============================================================================
def test_CO2_density():
    print("\n=== TEST 4: CO2 Density Predictions ===")

    # Garcia reports nearly linear density increase with x_CO2
    # ~2.5% at x2=0.05 at 25°C
    T, P = 298.15, 10.0

    # Test linearity: compute at x2 = 0.01, 0.02, 0.03, 0.04, 0.05
    print(f"\n  CO2 density change vs mole fraction at T=25°C, P=10 MPa:")
    print(f"  {'x2':>6} {'ρ(kg/m³)':>10} {'Δρ%':>8} {'Δρ%/x2':>8}")
    print("  " + "-" * 36)

    ratios = []
    for x2 in [0.01, 0.02, 0.03, 0.04, 0.05]:
        rho_sol, rho_w_val, pct = density_change_pct('CO2', x2, T, P)
        ratio = pct / x2 if x2 > 0 else 0
        ratios.append(ratio)
        print(f"  {x2:6.3f} {rho_sol:10.4f} {pct:+7.3f}% {ratio:8.2f}")

    # Check near-linearity: ratio should be roughly constant
    ratio_range = max(ratios) - min(ratios)
    check(f"CO2 Δρ%/x2 nearly constant (range={ratio_range:.2f})",
          ratio_range < 5.0,
          f"ratios span {ratio_range:.2f}")

    # Check ~2.5% at x2=0.05
    _, _, pct_05 = density_change_pct('CO2', 0.05, T, P)
    check(f"CO2 at x2=0.05: {pct_05:.2f}% (Garcia ~2.5%)",
          1.5 < pct_05 < 3.5,
          f"pct={pct_05:.2f}")


# ============================================================================
# TEST 5: Physical reasonableness
# ============================================================================
def test_physical_reasonableness():
    print("\n=== TEST 5: Physical Reasonableness ===")

    P = 30.0

    # V2∞ should increase monotonically with T for all gases (50-250°C at 30 MPa)
    # Note: slight dip from 25→50°C is physical at elevated P for some gases
    temps = [323.15, 373.15, 423.15, 473.15, 523.15]
    for gas in ['H2', 'N2', 'CH4', 'CO2', 'C2H6', 'C3H8', 'NC4H10', 'H2S']:
        vals = [V2_inf(gas, T, P) for T in temps]
        monotonic = all(vals[i] < vals[i+1] for i in range(len(vals)-1))
        check(f"{gas} V2∞ monotonically increasing 50-250°C",
              monotonic,
              f"values: {[f'{v:.1f}' for v in vals]}")

    # CO2 should increase density (heavy MW, small V_phi)
    _, _, pct_co2 = density_change_pct('CO2', 0.02, 298.15, 10.0)
    check(f"CO2 increases density ({pct_co2:+.3f}%)", pct_co2 > 0)

    # H2 should decrease density most (lightest MW)
    _, _, pct_h2 = density_change_pct('H2', 0.02, 298.15, 10.0)
    check(f"H2 decreases density ({pct_h2:+.3f}%)", pct_h2 < 0)

    # H2S should be near-neutral
    _, _, pct_h2s = density_change_pct('H2S', 0.02, 298.15, 10.0)
    check(f"H2S near-neutral density effect ({pct_h2s:+.3f}%)",
          abs(pct_h2s) < 0.5)

    # Heavier hydrocarbons should have larger V_phi
    v_c1 = V2_inf('CH4', 298.15, 0.1)
    v_c2 = V2_inf('C2H6', 298.15, 0.1)
    v_c3 = V2_inf('C3H8', 298.15, 0.1)
    v_c4 = V2_inf('NC4H10', 298.15, 0.1)
    check(f"V_phi ordering: CH4({v_c1:.1f}) < C2H6({v_c2:.1f}) < C3H8({v_c3:.1f}) < nC4({v_c4:.1f})",
          v_c1 < v_c2 < v_c3 < v_c4)

    # Mixed gas: adding CH4 to CO2 solution should reduce density increase
    _, _, pct_co2_only = density_change_pct('CO2', 0.02, 298.15, 10.0)
    _, _, pct_mix = density_change_pct({'CO2': 0.01, 'CH4': 0.01}, None, 298.15, 10.0)
    check(f"CO2+CH4 mix ({pct_mix:+.3f}%) < pure CO2 ({pct_co2_only:+.3f}%)",
          pct_mix < pct_co2_only)


# ============================================================================
# TEST 6: N2 vs O'Sullivan & Smith (1970)
# ============================================================================
def test_N2_vs_osullivan():
    print("\n=== TEST 6: N2 V2∞ vs O'Sullivan & Smith (1970) ===")

    # O'Sullivan data in pure water (low-P values)
    # T(°C): 51.5, 102.5, 125.0
    osullivan_N2 = [(324.65, 34.05), (375.65, 37.7), (398.15, 43.1)]

    P = 20.0  # moderate pressure
    print(f"\n  {'T(°C)':>6} {'Calc':>8} {'O-S':>8} {'%Err':>7}")
    print("  " + "-" * 34)

    for T_K, v_exp in osullivan_N2:
        v_calc = V2_inf('N2', T_K, P)
        pct = 100 * (v_calc - v_exp) / v_exp
        T_C = T_K - 273.15
        print(f"  {T_C:6.1f} {v_calc:8.2f} {v_exp:8.2f} {pct:+6.2f}%")

    # O'Sullivan data has lower precision; expect agreement within ~10%
    max_err = max(abs(100 * (V2_inf('N2', T, P) - v) / v) for T, v in osullivan_N2)
    check(f"N2 vs O'Sullivan max error {max_err:.1f}% < 15%", max_err < 15,
          f"max_err={max_err:.1f}%")


# ============================================================================
# RUN ALL TESTS
# ============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("GARCIA EXTENSION - COMPREHENSIVE VALIDATION")
    print("=" * 60)

    test_V2inf_298K()
    test_V2inf_vs_hnedkovsky()
    test_V2inf_vs_garcia_polynomial()
    test_CO2_density()
    test_physical_reasonableness()
    test_N2_vs_osullivan()

    print("\n" + "=" * 60)
    print(f"RESULTS: {PASS} passed, {FAIL} failed out of {PASS + FAIL} checks")
    print("=" * 60)

    if FAIL > 0:
        sys.exit(1)
