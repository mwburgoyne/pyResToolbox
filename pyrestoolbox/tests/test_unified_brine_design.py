#!/usr/bin/env python3
"""
Test scripts for the unified brine density & viscosity framework design.

These are exploratory/validation tests, NOT part of the main test suite yet.
They help understand the current implementations and validate design decisions.

Run with: PYTHONPATH=/home/mark/projects python3 tests/test_unified_brine_design.py
"""

import sys
import numpy as np

sys.path.insert(0, '/home/mark/projects')

from pyrestoolbox import brine

# ============================================================================
# R12: Compare Spivey Pitzer CH4 density vs Garcia-form CH4 density
# ============================================================================
# The Spivey/McCain approach (brine_props) computes CH4-saturated brine density
# using Pitzer interaction parameter derivatives (Eq 4.22-4.24).
#
# The Garcia approach uses V_phi(T) polynomial + Eq 18 mixing rule.
# For CO2: V_phi = 37.51 - 0.09585*T + 8.74e-4*T^2 - 5.044e-7*T^3
#
# We need to understand: what effective V_phi does the Spivey method imply
# for CH4, and how does it compare to literature values?
# ============================================================================

def extract_spivey_ch4_vphi():
    """
    Extract the effective apparent molar volume of dissolved CH4
    implied by the Spivey/McCain approach at various temperatures.

    The Spivey method computes:
        Vmch4b = R*T * (du/dp + 2*m*dlambda/dp + m^2*deta/dp)  [Eq 4.22]
    where Vmch4b is partial molar volume in cm3/mol at the given T,P,salinity.

    We'll extract Vmch4b at various conditions and compare to published
    V_phi^inf values for CH4 in water.
    """
    from pyrestoolbox.brine.brine import _Eq41, _RHOW_T70_ARR, _EWT_ARR, _FWT_ARR
    from pyrestoolbox.brine.brine import _DM2T_ARR, _DM32T_ARR, _DM1T_ARR, _DM12T_ARR
    from pyrestoolbox.brine.brine import _EMT_ARR, _FM32T_ARR, _FM1T_ARR, _FM12T_ARR

    print("=" * 80)
    print("R12: Spivey Pitzer CH4 Partial Molar Volume vs Temperature")
    print("=" * 80)
    print()

    # Spivey Pitzer interaction parameter coefficients for CH4
    u_arr = [0, 8.3143711, -7.2772168e-4, 2.1489858e3, -1.4019672e-5,
             -6.6743449e5, 7.698589e-2, -5.0253331e-5, -30.092013, 4.8468502e3, 0]
    lambda_arr = [0, -0.80898, 1.0827e-3, 183.85, 0, 0, 3.924e-4, 0, 0, 0, -1.97e-6]

    print(f"{'T(°F)':>8} {'T(°C)':>8} {'P(psia)':>8} {'Wt%':>6} {'Vmch4b':>10} {'rho_brine':>10} {'rho_ch4sat':>12} {'delta_rho%':>10}")
    print("-" * 90)

    results = []

    for degf in [100, 150, 200, 250, 300]:
        for p in [2000, 5000, 8000]:
            for wt in [0, 3, 10]:
                degc = (degf - 32) / 1.8
                degk = degc + 273
                Mpa = p * 0.00689476
                m = 1000 * (wt / 100) / (58.4428 * (1 - wt / 100)) if wt > 0 else 0

                # Compute du/dp and dlambda/dp (pressure derivatives of Pitzer params)
                dudptm = (u_arr[6] + u_arr[7] * degk + u_arr[8] / degk + u_arr[9] / (degk * degk))
                dlambdadptm = lambda_arr[6] + 2 * lambda_arr[10] * Mpa

                # Vmch4b = R*T*(du/dp + 2*m*dlambda/dp)  [cm3/mol]
                Vmch4b = 8.314467 * degk * (dudptm + 2 * m * dlambdadptm)

                # Also get full brine_props results for density comparison
                bw, lden, visw, cw, rsw = brine.brine_props(p=p, degf=degf, wt=wt, ch4_sat=1.0)
                bw0, lden0, visw0, cw0, rsw0 = brine.brine_props(p=p, degf=degf, wt=wt, ch4_sat=0.0)

                delta_pct = (lden - lden0) / lden0 * 100

                results.append({
                    'degf': degf, 'degc': degc, 'p': p, 'wt': wt,
                    'Vmch4b': Vmch4b, 'rho_nogas': lden0, 'rho_ch4sat': lden,
                    'delta_pct': delta_pct
                })

                if wt == 0 and p == 5000:  # Print subset
                    print(f"{degf:>8.0f} {degc:>8.1f} {p:>8.0f} {wt:>6.0f} {Vmch4b:>10.2f} {lden0:>10.5f} {lden:>12.5f} {delta_pct:>10.3f}")

    print()
    print("Key observations:")
    print("  - Vmch4b is the Spivey/Pitzer partial molar volume of CH4 in brine (cm3/mol)")
    print("  - Garcia CO2 V_phi at 25°C = 37.51 cm3/mol; literature CH4 V_phi ~37 cm3/mol")
    print("  - delta_rho% shows density change from dissolved CH4 (should be NEGATIVE for CH4)")
    print()

    # Summary statistics
    vmch4_values = [r['Vmch4b'] for r in results if r['wt'] == 0]
    print(f"  Vmch4b range (fresh water): {min(vmch4_values):.2f} to {max(vmch4_values):.2f} cm3/mol")

    delta_values = [r['delta_pct'] for r in results]
    print(f"  Density change range: {min(delta_values):.3f}% to {max(delta_values):.3f}%")

    return results


def compare_garcia_co2_implementations():
    """
    Verify that the Garcia CO2 density correction in CO2_Brine_Mixture
    produces consistent results. This is validation for TASK D6.
    """
    print()
    print("=" * 80)
    print("D6: Garcia CO2 density implementation verification")
    print("=" * 80)
    print()

    print(f"{'T(°C)':>8} {'P(bar)':>8} {'ppm':>8} {'xCO2':>10} {'rho_co2sat':>12} {'rho_brine':>10} {'delta%':>8}")
    print("-" * 75)

    test_cases = [
        (100, 85, 0),      # Moderate conditions, fresh water
        (200, 85, 0),      # Higher pressure
        (300, 85, 0),      # High pressure
        (100, 85, 30000),  # 3% NaCl
        (200, 85, 30000),
        (300, 150, 0),     # High T high P
        (200, 150, 50000), # High salinity
    ]

    for pbar, degc, ppm in test_cases:
        try:
            mix = brine.CO2_Brine_Mixture(pres=pbar, temp=degc, ppm=ppm, metric=True)
            xCO2 = mix.x[0]
            rho_sat = mix.bDen[0]
            rho_brine = mix.bDen[1]
            delta = (rho_sat - rho_brine) / rho_brine * 100

            print(f"{degc:>8.0f} {pbar:>8.0f} {ppm:>8.0f} {xCO2:>10.5f} {rho_sat:>12.5f} {rho_brine:>10.5f} {delta:>8.3f}")
        except Exception as e:
            print(f"{degc:>8.0f} {pbar:>8.0f} {ppm:>8.0f} {'ERROR':>10} {str(e)[:40]}")

    print()
    print("Key observations:")
    print("  - CO2 INCREASES brine density (positive delta)")
    print("  - Effect is ~1-3% at typical reservoir conditions")
    print("  - Effect decreases with temperature (V_phi increases)")


def compare_viscosity_corrections():
    """
    Compare viscosity with and without dissolved gas corrections.
    Shows magnitude of effect and validates the Islam-Carlson CO2 correction.
    """
    print()
    print("=" * 80)
    print("Viscosity correction magnitude comparison")
    print("=" * 80)
    print()

    print("--- CO2 effect on brine viscosity (Islam-Carlson) ---")
    print(f"{'T(°C)':>8} {'P(bar)':>8} {'ppm':>8} {'xCO2':>10} {'vis_co2':>10} {'vis_brine':>10} {'increase%':>10}")
    print("-" * 75)

    for pbar, degc, ppm in [(200, 85, 0), (200, 85, 30000), (300, 85, 0), (300, 150, 0), (100, 50, 0)]:
        try:
            mix = brine.CO2_Brine_Mixture(pres=pbar, temp=degc, ppm=ppm, metric=True)
            xCO2 = mix.x[0]
            vis_sat = mix.bVis[0]
            vis_brine = mix.bVis[1]
            increase = (vis_sat - vis_brine) / vis_brine * 100

            print(f"{degc:>8.0f} {pbar:>8.0f} {ppm:>8.0f} {xCO2:>10.5f} {vis_sat:>10.4f} {vis_brine:>10.4f} {increase:>10.2f}")
        except Exception as e:
            print(f"{degc:>8.0f} {pbar:>8.0f} {ppm:>8.0f} {'ERROR':>10} {str(e)[:40]}")

    print()
    print("--- CH4 effect on brine viscosity (currently: NONE in Spivey) ---")
    print(f"{'T(°F)':>8} {'P(psi)':>8} {'wt%':>6} {'vis_ch4sat':>12} {'vis_nogas':>10} {'change%':>10}")
    print("-" * 60)

    for degf, p, wt in [(200, 5000, 0), (200, 5000, 3), (250, 5000, 0), (150, 3000, 0)]:
        bw1, den1, vis1, cw1, rsw1 = brine.brine_props(p=p, degf=degf, wt=wt, ch4_sat=1.0)
        bw0, den0, vis0, cw0, rsw0 = brine.brine_props(p=p, degf=degf, wt=wt, ch4_sat=0.0)
        change = (vis1 - vis0) / vis0 * 100

        print(f"{degf:>8.0f} {p:>8.0f} {wt:>6.0f} {vis1:>12.4f} {vis0:>10.4f} {change:>10.2f}")

    print()
    print("Key observations:")
    print("  - CO2: Islam-Carlson gives ~5-15% viscosity increase at saturation")
    print("  - CH4: Currently 0% change (Spivey doesn't correct for dissolved CH4)")
    print("  - Literature suggests CH4 effect is ~6% max (SPE-14211)")
    print("  - For N2, H2: effect is negligible due to low solubility")


def demonstrate_garcia_generalization():
    """
    Demonstrate how the Garcia mixing rule generalizes to any gas.
    Uses the existing CO2 implementation as the template and shows
    what the unified function signature would look like.
    """
    print()
    print("=" * 80)
    print("Demonstration: Garcia mixing rule with different gases")
    print("=" * 80)
    print()

    # Garcia Eq 18: rho = (1 + MW_gas/MW_brine * x/(1-x)) / (V_phi * x/((1-x)*MW_brine) + 1/rho_brine)
    # This is IDENTICAL for any gas — only V_phi and MW change

    def garcia_density(rho_brine, T_celsius, x_gas, MW_brine, MW_gas, vphi_fn):
        """Generic Garcia Eq 18"""
        x_not_gas = 1.0 - x_gas
        x_ratio = x_gas / x_not_gas
        mw_ratio = MW_gas / MW_brine
        vphi = vphi_fn(T_celsius)
        return (1.0 + mw_ratio * x_ratio) / (vphi * x_ratio / MW_brine + 1.0 / rho_brine)

    # V_phi polynomials: V_phi(T) = a0 + a1*T + a2*T^2 + a3*T^3 where T in °C
    # CO2: from Garcia (2001) — these are the KNOWN coefficients
    vphi_co2 = lambda T: 37.51 + T * (-0.09585 + T * (8.74e-4 + T * (-5.044e-7)))

    # CH4: PLACEHOLDER — approximate from literature ~37 cm3/mol at 25°C
    # Real coefficients need to be fitted from Hnedkovsky (1996) or Plyasunov (2019)
    # For now, use a rough estimate based on similarity to CO2 V_phi shape
    vphi_ch4_placeholder = lambda T: 37.0 + T * (0.10) + T**2 * (1.0e-4)  # PLACEHOLDER ONLY

    # H2S: PLACEHOLDER — V_phi ~35 cm3/mol at 25°C (from extending_garcia doc)
    vphi_h2s_placeholder = lambda T: 35.0 + T * (-0.05)  # PLACEHOLDER ONLY

    MW = {'CO2': 44.01, 'CH4': 16.043, 'H2S': 34.082, 'N2': 28.014, 'H2': 2.016}

    # Test conditions
    T_celsius = 85.0
    rho_brine = 1.005  # typical brine density g/cm3
    MW_brine = 18.15   # slightly saline brine MW

    print(f"Conditions: T = {T_celsius}°C, rho_brine = {rho_brine} g/cm3")
    print()
    print(f"{'Gas':>6} {'MW':>8} {'V_phi':>8} {'MW/V_phi':>8} {'x_gas':>8} {'rho_sat':>10} {'delta%':>8} {'Effect':>12}")
    print("-" * 80)

    for gas_name, vphi_fn, mw in [
        ('CO2', vphi_co2, MW['CO2']),
        ('CH4', vphi_ch4_placeholder, MW['CH4']),
        ('H2S', vphi_h2s_placeholder, MW['H2S']),
    ]:
        x_gas = 0.02  # typical dissolved mole fraction
        vphi_val = vphi_fn(T_celsius)
        mw_over_vphi = mw / vphi_val
        rho_sat = garcia_density(rho_brine, T_celsius, x_gas, MW_brine, mw, vphi_fn)
        delta = (rho_sat - rho_brine) / rho_brine * 100
        effect = "increases" if delta > 0 else "decreases"

        print(f"{gas_name:>6} {mw:>8.3f} {vphi_val:>8.2f} {mw_over_vphi:>8.3f} {x_gas:>8.4f} {rho_sat:>10.5f} {delta:>8.3f} {effect:>12}")

    print()
    print("Note: CH4 and H2S V_phi values are PLACEHOLDERS. Real coefficients needed.")
    print("Note: MW/V_phi > ~1.0 → increases density; MW/V_phi < ~1.0 → decreases density")
    print()

    # Mixed gas demonstration
    print("--- Mixed gas example ---")
    print("If multiple gases dissolved, use mole-fraction-weighted V_phi and MW:")

    x_co2 = 0.015
    x_ch4 = 0.005
    x_total = x_co2 + x_ch4
    y_co2 = x_co2 / x_total  # fraction of CO2 among dissolved gases
    y_ch4 = x_ch4 / x_total

    vphi_eff = y_co2 * vphi_co2(T_celsius) + y_ch4 * vphi_ch4_placeholder(T_celsius)
    mw_eff = y_co2 * MW['CO2'] + y_ch4 * MW['CH4']

    vphi_eff_fn = lambda T: vphi_eff  # constant for this demo
    rho_mixed = garcia_density(rho_brine, T_celsius, x_total, MW_brine, mw_eff, vphi_eff_fn)
    delta_mixed = (rho_mixed - rho_brine) / rho_brine * 100

    print(f"  CO2: x={x_co2}, CH4: x={x_ch4}, total: x={x_total}")
    print(f"  V_phi_eff = {vphi_eff:.2f} cm3/mol, MW_eff = {mw_eff:.3f} g/mol")
    print(f"  rho_mixed = {rho_mixed:.5f} g/cm3 (delta = {delta_mixed:.3f}%)")


def print_unified_api_design():
    """Print the proposed unified API design."""
    print()
    print("=" * 80)
    print("PROPOSED UNIFIED API DESIGN")
    print("=" * 80)
    print("""
Option A: Single function approach
-----------------------------------
brine.brine_density(
    p,              # Pressure (psia)
    degf,           # Temperature (°F)
    wt=0,           # Salt wt% (0-100)
    x_ch4=0,        # Dissolved CH4 mole fraction in brine
    x_co2=0,        # Dissolved CO2 mole fraction in brine
    x_h2s=0,        # Dissolved H2S mole fraction in brine
    x_n2=0,         # Dissolved N2 mole fraction in brine
    x_h2=0,         # Dissolved H2 mole fraction in brine
    x_c2h6=0,       # Dissolved C2H6 mole fraction in brine
    x_c3h8=0,       # Dissolved C3H8 mole fraction in brine
) -> float          # Density in g/cm3

brine.brine_viscosity(
    p, degf, wt=0,
    x_ch4=0, x_co2=0, x_h2s=0, x_n2=0, x_h2=0, x_c2h6=0, x_c3h8=0,
) -> float          # Viscosity in cP

brine.brine_properties(
    p, degf, wt=0,
    x_ch4=0, x_co2=0, x_h2s=0, x_n2=0, x_h2=0, x_c2h6=0, x_c3h8=0,
) -> dict           # All properties: density, viscosity, Bw, compressibility

Notes:
- Gas-free brine density from Spivey (shared engine, already exists)
- Dissolved gas density correction from Garcia Eq 18 with per-gas V_phi(T)
- Base viscosity from Mao-Duan (already exists)
- Viscosity correction: multiplicative per-gas (Islam-Carlson form)
- Dissolved gas mole fractions come from the CALLER (could be from
  Spycher-Pruess for CO2, Spivey for CH4, future Soreide-Whitson for
  multi-gas, or any external source)
- This DECOUPLES density/viscosity from solubility models

Option B: Keep existing APIs, add unified function alongside
------------------------------------------------------------
- Keep brine_props() as-is (backward compatible, CH4-focused)
- Keep CO2_Brine_Mixture as-is (backward compatible, CO2-focused)
- Add new unified functions that accept arbitrary compositions
- Existing functions could internally call the unified functions
  once validated

RECOMMENDED: Option B for initial implementation
  - No breaking changes
  - Allows side-by-side validation
  - Migrate users gradually
  - Eventually deprecate old functions
""")


if __name__ == '__main__':
    # Run all explorations
    extract_spivey_ch4_vphi()
    compare_garcia_co2_implementations()
    compare_viscosity_corrections()
    demonstrate_garcia_generalization()
    print_unified_api_design()
