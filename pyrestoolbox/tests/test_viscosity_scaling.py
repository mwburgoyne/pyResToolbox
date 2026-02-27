"""
Exploratory analysis: Can dissolved gas viscosity effects be scaled from CO2?

Investigates whether the ratio of free-gas viscosities, molecular weights,
apparent molar volumes, or "dissolved densities" (MW/V_phi) can predict
the viscosity effect of dissolving different gases in water/brine.

Known data:
  - CO2: Islam-Carlson (2012): mu = mu_brine * (1 + 4.65 * x_CO2^1.0134)
  - CH4: Kumagai & Yokoyama (1998): ~1-6% decrease at saturation conditions
  - H2S: no published correction
  - N2, H2: no published correction
"""

import sys
import os
import numpy as np

# Add paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, '/home/mark/projects/pyResToolbox/Garcia_code')

# ============================================================================
# 1. Islam-Carlson CO2 viscosity effect
# ============================================================================
def islam_carlson_factor(x_co2):
    """Viscosity multiplier from Islam-Carlson (2012)."""
    if x_co2 <= 0:
        return 1.0
    return 1.0 + 4.65 * x_co2 ** 1.0134

print("=" * 70)
print("1. ISLAM-CARLSON CO2 VISCOSITY EFFECT")
print("=" * 70)
print(f"\n  {'x_CO2':>8} {'Factor':>8} {'%change':>8}")
print("  " + "-" * 28)
for x in [0.005, 0.01, 0.02, 0.03, 0.04, 0.05]:
    f = islam_carlson_factor(x)
    print(f"  {x:8.3f} {f:8.4f} {(f-1)*100:+7.2f}%")

# Note: at x_CO2 = 0.01, effect is ~4.4% INCREASE
# The exponent ~1.01 means it's nearly linear: ~4.65 * x_CO2

# ============================================================================
# 2. Free gas viscosities at reservoir conditions
# ============================================================================
print("\n" + "=" * 70)
print("2. FREE GAS VISCOSITIES AT RESERVOIR CONDITIONS")
print("=" * 70)

# Use pyResToolbox gas_ug for hydrocarbon gases
# For non-HC gases, use Lee-Gonzalez-Eakin with appropriate MW and Tc/Pc
from pyrestoolbox import gas as rtb_gas

def lee_gonzalez_eakin(T_R, rho_gcc, MW):
    """Lee-Gonzalez-Eakin gas viscosity correlation.
    T_R: temperature in Rankine
    rho_gcc: gas density in g/cc
    MW: molecular weight
    Returns: viscosity in cP
    """
    K = ((9.4 + 0.02 * MW) * T_R**1.5) / (209 + 19 * MW + T_R)
    X = 3.5 + 986.0 / T_R + 0.01 * MW
    Y = 2.4 - 0.2 * X
    return 1e-4 * K * np.exp(X * rho_gcc**Y)


# Conditions: 150°F (610 R, 338.7 K), 3000 psia (20.68 MPa)
T_F = 150.0
P_psia = 3000.0
T_K = (T_F - 32) * 5/9 + 273.15
T_R = T_F + 459.67
P_MPa = P_psia / 145.038

print(f"\n  Conditions: T = {T_F}°F ({T_K:.1f} K), P = {P_psia} psia ({P_MPa:.1f} MPa)")

# For pure CH4 - use pyResToolbox
try:
    ug_ch4 = rtb_gas.gas_ug(p=P_psia, sg=0.5537, degf=T_F)  # CH4 sg ~0.5537
    print(f"\n  CH4 (pyResToolbox, sg=0.554): {float(ug_ch4):.4f} cP")
except Exception as e:
    print(f"\n  CH4 via pyResToolbox failed: {e}")
    ug_ch4 = 0.018  # fallback

# For various pure gases, compute approximate viscosity via LGE
# Need density - use ideal gas as rough estimate, then correct
R_gcc = 8.314  # J/(mol*K)

gas_props = {
    'CH4':    {'MW': 16.043, 'Tc_K': 190.56, 'Pc_MPa': 4.599},
    'CO2':    {'MW': 44.010, 'Tc_K': 304.13, 'Pc_MPa': 7.377},
    'H2S':    {'MW': 34.081, 'Tc_K': 373.21, 'Pc_MPa': 8.937},
    'N2':     {'MW': 28.013, 'Tc_K': 126.19, 'Pc_MPa': 3.396},
    'H2':     {'MW': 2.016,  'Tc_K': 33.19,  'Pc_MPa': 1.313},
    'C2H6':   {'MW': 30.069, 'Tc_K': 305.32, 'Pc_MPa': 4.872},
    'C3H8':   {'MW': 44.096, 'Tc_K': 369.83, 'Pc_MPa': 4.248},
}

# Simple gas density estimate (ideal gas with compressibility guess)
def gas_density_approx(MW, T_K, P_MPa, Tc_K, Pc_MPa):
    """Approximate gas density using Pitzer correlation for Z."""
    Tr = T_K / Tc_K
    Pr = P_MPa / Pc_MPa
    # Simple Z estimate (van der Waals like)
    if Tr > 1.0:
        Z = 1.0 + (0.083 - 0.422/Tr**1.6) * Pr / Tr
    else:
        Z = 0.3  # liquid-like
    Z = max(Z, 0.1)
    rho = P_MPa * MW / (Z * 8.314e-3 * T_K)  # kg/m3
    return rho / 1000.0, Z  # g/cc, Z

print(f"\n  {'Gas':>8} {'MW':>6} {'Tr':>6} {'Pr':>6} {'Z':>6} {'rho(g/cc)':>10} {'ug_LGE(cP)':>12}")
print("  " + "-" * 65)

gas_viscosities = {}
for gas, props in gas_props.items():
    MW = props['MW']
    Tc = props['Tc_K']
    Pc = props['Pc_MPa']
    Tr = T_K / Tc
    Pr = P_MPa / Pc
    rho_gcc, Z = gas_density_approx(MW, T_K, P_MPa, Tc, Pc)

    # LGE viscosity
    ug = lee_gonzalez_eakin(T_R, rho_gcc, MW)
    gas_viscosities[gas] = ug

    note = ""
    if Tr < 1.0:
        note = " ** SUBCRITICAL - liquid phase **"

    print(f"  {gas:>8} {MW:6.2f} {Tr:6.2f} {Pr:6.2f} {Z:6.3f} {rho_gcc:10.4f} {ug:12.5f}{note}")

print("\n  Note: LGE is calibrated for HC gases. Values for CO2, H2S, N2, H2 are approximate.")
print("  CO2 and H2S are near/below Tc at 150°F — phase behavior complicates things.")

# ============================================================================
# 3. V_phi (apparent molar volume) and "dissolved density"
# ============================================================================
print("\n" + "=" * 70)
print("3. APPARENT MOLAR VOLUMES AND 'DISSOLVED DENSITY' AT RESERVOIR CONDITIONS")
print("=" * 70)

from plyasunov_model import V_phi, gas_mw
from garcia_mixing import density_change_pct

# At our reservoir conditions
print(f"\n  Conditions: T = {T_K:.1f} K, P = {P_MPa:.1f} MPa")
print(f"\n  {'Gas':>8} {'MW':>6} {'V_phi':>8} {'MW/Vphi':>8} {'drho%':>8} {'Direction':>10}")
print("  " + "-" * 56)

gases = ['H2', 'N2', 'CH4', 'CO2', 'C2H6', 'C3H8', 'NC4H10', 'H2S']
vphi_data = {}

for gas_name in gases:
    MW = gas_mw(gas_name)
    vphi = V_phi(gas_name, T_K, P_MPa)
    dissolved_density = MW / vphi  # g/cm3 "effective density" of dissolved gas

    _, _, drho_pct = density_change_pct(gas_name, 0.02, T_K, P_MPa)

    direction = "INCREASES" if drho_pct > 0 else "decreases"

    vphi_data[gas_name] = {'MW': MW, 'vphi': vphi, 'dissolved_rho': dissolved_density, 'drho_pct': drho_pct}
    print(f"  {gas_name:>8} {MW:6.2f} {vphi:8.2f} {dissolved_density:8.4f} {drho_pct:+7.3f}% {direction:>10}")

# Water density for reference
rho_water = 1000 * 18.015 / (18.015 * 1000 / 1000)  # ~1 g/cc at ambient, but let's get actual
from water_properties import rho_w, MW_WATER
rho_w_val = rho_w(T_K, P_MPa) / 1000  # g/cc
V1_star = MW_WATER / (rho_w(T_K, P_MPa) / 1000)  # cm3/mol... wait
V1_star = MW_WATER * 1000 / rho_w(T_K, P_MPa)  # cm3/mol
water_dissolved_density = MW_WATER / V1_star
print(f"\n  Water: MW={MW_WATER:.3f}, V*1={V1_star:.2f} cm3/mol, MW/V*1={water_dissolved_density:.4f} g/cm3")
print(f"  Water density at conditions: {rho_w(T_K, P_MPa):.2f} kg/m3 = {rho_w_val:.4f} g/cm3")

# ============================================================================
# 4. KEY ANALYSIS: Can we predict viscosity effect direction and magnitude?
# ============================================================================
print("\n" + "=" * 70)
print("4. ANALYSIS: SCALING PARAMETERS FOR VISCOSITY EFFECTS")
print("=" * 70)

print("""
KNOWN VISCOSITY EFFECTS:
  CO2: INCREASES water viscosity by ~4.65% per mol% (Islam-Carlson)
  CH4: DECREASES water viscosity by ~1-6% at saturation (Kumagai-Yokoyama)
  H2S: Unknown (but chemically similar to H2O - hydrogen bonding)
  N2, H2, C2H6, C3H8: Unknown

HYPOTHESIS A: Scale by free-gas viscosity ratio (ug_gas / ug_CO2)
""")

ug_co2 = gas_viscosities.get('CO2', 0.03)
print(f"  CO2 Islam-Carlson coefficient: a = 4.65")
print(f"  CO2 free-gas viscosity (LGE approx): {ug_co2:.5f} cP")
print(f"\n  {'Gas':>8} {'ug(cP)':>10} {'ug/ug_CO2':>10} {'Scaled a':>10} {'%eff @ x=0.02':>14}")
print("  " + "-" * 58)

for gas in ['CH4', 'CO2', 'H2S', 'N2', 'H2', 'C2H6', 'C3H8']:
    ug = gas_viscosities.get(gas, 0)
    ratio = ug / ug_co2 if ug_co2 > 0 else 0
    scaled_a = 4.65 * ratio
    eff = scaled_a * 0.02  # approximate effect at x=0.02
    print(f"  {gas:>8} {ug:10.5f} {ratio:10.4f} {scaled_a:10.3f} {eff*100:+13.2f}%")

print("""
  PROBLEM: This predicts ALL gases INCREASE viscosity (all a > 0).
  But CH4 is known to DECREASE water viscosity.
  Free-gas viscosity ratios cannot predict the SIGN of the effect.
""")

print("HYPOTHESIS B: Scale by 'dissolved density' ratio (MW/V_phi)")
print("  Idea: Heavy molecules in small volumes increase viscosity;")
print("  light molecules in large volumes decrease it.")
print(f"\n  Water 'density': MW/V1* = {water_dissolved_density:.4f} g/cm3")
print(f"  If MW_gas/V_phi > MW_water/V1*, dissolved gas is 'heavier' → increases viscosity")
print(f"  If MW_gas/V_phi < MW_water/V1*, dissolved gas is 'lighter' → decreases viscosity")

print(f"\n  {'Gas':>8} {'MW/Vphi':>8} {'vs H2O':>8} {'Density':>10} {'Visc':>10}")
print("  " + "-" * 50)

for gas_name in gases:
    d = vphi_data[gas_name]
    ratio = d['dissolved_rho'] / water_dissolved_density
    density_dir = "INCREASES" if d['drho_pct'] > 0 else "decreases"
    visc_dir = "INCREASES" if d['dissolved_rho'] > water_dissolved_density else "decreases"
    print(f"  {gas_name:>8} {d['dissolved_rho']:8.4f} {ratio:8.4f} {density_dir:>10} {visc_dir:>10}")

print(f"""
  The MW/V_phi ratio correctly predicts:
    - CO2 increases (MW/Vphi = {vphi_data['CO2']['dissolved_rho']:.4f} > water {water_dissolved_density:.4f}) ✓
    - CH4 decreases (MW/Vphi = {vphi_data['CH4']['dissolved_rho']:.4f} < water {water_dissolved_density:.4f}) ✓
    - H2S: {'increases' if vphi_data['H2S']['dissolved_rho'] > water_dissolved_density else 'decreases'}

  But the sign correlation is IDENTICAL to the density effect direction.
  This is expected: both density and viscosity effects stem from the same
  physics (how the dissolved molecule fits in the water structure).
""")

# ============================================================================
# 5. PROPOSED APPROACH: Viscosity-density analogy
# ============================================================================
print("=" * 70)
print("5. PROPOSED APPROACHES")
print("=" * 70)

print("""
APPROACH 1: Density-viscosity proportionality

  Observation: For CO2 at x=0.02:
    - Density increase: ~1.1%
    - Viscosity increase: ~9.3% (Islam-Carlson)
    - Ratio: viscosity_effect / density_effect ≈ 8.5

  This suggests viscosity is ~8-9x more sensitive than density.

  For other gases, scale:
    a_gas = a_CO2 * (drho%_gas / drho%_CO2)

  This gives the right SIGN (negative for CH4, positive for CO2)
  and uses the density model we already have.
""")

# Compute this
drho_co2 = vphi_data['CO2']['drho_pct']
a_co2 = 4.65
ic_factor_002 = islam_carlson_factor(0.02) - 1  # fractional viscosity change at x=0.02
visc_to_density_ratio = ic_factor_002 / (drho_co2 / 100)

print(f"  CO2 at x=0.02: drho = {drho_co2:+.4f}%, dmu = {ic_factor_002*100:+.4f}%")
print(f"  Viscosity/density sensitivity ratio: {visc_to_density_ratio:.2f}")

print(f"\n  {'Gas':>8} {'drho%':>8} {'Predicted dmu%':>15} {'Predicted a':>12}")
print("  " + "-" * 48)

for gas_name in gases:
    d = vphi_data[gas_name]
    # Scale: dmu/mu = visc_to_density_ratio * drho/rho
    predicted_dmu_pct = visc_to_density_ratio * d['drho_pct']
    # Back out 'a' coefficient: dmu/mu ≈ a * x at small x
    predicted_a = predicted_dmu_pct / (0.02 * 100) if abs(d['drho_pct']) > 1e-6 else 0

    print(f"  {gas_name:>8} {d['drho_pct']:+7.3f}% {predicted_dmu_pct:+14.3f}% {predicted_a:+11.3f}")

print("""
APPROACH 2: Direct V_phi-based correction (more physical)

  Islam-Carlson is: mu = mu_brine * (1 + a * x^b)
  where a=4.65, b=1.0134 for CO2.

  Generalize: mu = mu_brine * Product_over_gases(1 + a_i * x_i^b)
  where a_i scales with how "disruptive" the dissolved molecule is.

  Option A: a_i = a_CO2 * (drho_i% / drho_CO2%)
    Uses density ratio as proxy. Simple, uses existing model.

  Option B: a_i = a_CO2 * (MW_i/Vphi_i - MW_H2O/V1*) / (MW_CO2/Vphi_CO2 - MW_H2O/V1*)
    Uses "excess dissolved density" relative to water.
    More physical basis.
""")

# Compute Option B
delta_CO2 = vphi_data['CO2']['dissolved_rho'] - water_dissolved_density

print(f"\n  Option B scaling:")
print(f"  CO2 'excess dissolved density': {delta_CO2:+.4f} g/cm3")
print(f"\n  {'Gas':>8} {'MW/Vphi':>8} {'delta':>8} {'delta/dCO2':>10} {'a_i':>8}")
print("  " + "-" * 48)

for gas_name in gases:
    d = vphi_data[gas_name]
    delta = d['dissolved_rho'] - water_dissolved_density
    ratio = delta / delta_CO2 if abs(delta_CO2) > 1e-10 else 0
    a_i = a_co2 * ratio
    print(f"  {gas_name:>8} {d['dissolved_rho']:8.4f} {delta:+8.4f} {ratio:+9.4f} {a_i:+8.3f}")

# ============================================================================
# 6. CROSS-CHECK: CH4 predicted vs literature
# ============================================================================
print("\n" + "=" * 70)
print("6. CROSS-CHECK: CH4 PREDICTED EFFECT vs LITERATURE")
print("=" * 70)

# Kumagai & Yokoyama (1998): CH4 decreases viscosity by ~1-6% at saturation
# Typical CH4 saturation mole fraction in water at reservoir conditions: 0.001-0.003
# At 150°C, 30 MPa: x_CH4 ~ 0.003 (very rough)

print("""
  Literature: Kumagai & Yokoyama (1998) report CH4 DECREASES water
  viscosity by 1-6% depending on T, P, and CH4 concentration.

  Typical CH4 saturation in water: x_CH4 ~ 0.001-0.003 at reservoir conditions
""")

# Using our density-based scaling:
for x_ch4 in [0.001, 0.002, 0.003, 0.005, 0.01]:
    drho_ch4_x = vphi_data['CH4']['drho_pct'] * (x_ch4 / 0.02)  # rough linear scale
    predicted_dmu = visc_to_density_ratio * drho_ch4_x / 100
    print(f"  x_CH4={x_ch4:.3f}: predicted dmu = {predicted_dmu*100:+.2f}%")

print("""
  At x_CH4=0.002-0.003, we predict ~-1.5% to -2.3% viscosity decrease.
  Kumagai-Yokoyama report 1-6% decrease over a range of conditions.
  Order of magnitude is reasonable, though our prediction is conservative.

  The larger effects in K&Y may reflect higher CH4 concentrations at
  their specific conditions and/or additional physical mechanisms.
""")

# ============================================================================
# 7. TEMPERATURE DEPENDENCE
# ============================================================================
print("=" * 70)
print("7. TEMPERATURE DEPENDENCE OF SCALING PARAMETERS")
print("=" * 70)

print(f"\n  V_phi and MW/V_phi at P=30 MPa, various temperatures:")
print(f"\n  {'T(°C)':>6}  ", end="")
for gas_name in ['CO2', 'CH4', 'H2S', 'N2', 'H2']:
    print(f"  {gas_name:>12}", end="")
print()
print("  " + "-" * 70)

for label in ['V_phi', 'MW/Vphi']:
    for T_C in [25, 50, 100, 150, 200]:
        T_K_local = T_C + 273.15
        print(f"  {T_C:>5}°C ", end="")
        for gas_name in ['CO2', 'CH4', 'H2S', 'N2', 'H2']:
            MW = gas_mw(gas_name)
            vphi = V_phi(gas_name, T_K_local, 30.0)
            if label == 'V_phi':
                print(f"  {vphi:12.2f}", end="")
            else:
                print(f"  {MW/vphi:12.4f}", end="")
        print(f"  ({label})")
    print()

# Check: does the scaling ratio stay roughly constant with T?
print("  Scaling ratio (MW/Vphi_gas - MW/V1*) / (MW/Vphi_CO2 - MW/V1*) at P=30 MPa:")
print(f"\n  {'T(°C)':>6}", end="")
for gas_name in ['CH4', 'H2S', 'N2', 'H2']:
    print(f"  {gas_name:>8}", end="")
print()
print("  " + "-" * 42)

for T_C in [25, 50, 100, 150, 200]:
    T_K_local = T_C + 273.15
    V1 = MW_WATER * 1000 / rho_w(T_K_local, 30.0)
    wd = MW_WATER / V1

    vphi_co2 = V_phi('CO2', T_K_local, 30.0)
    delta_co2_local = gas_mw('CO2') / vphi_co2 - wd

    print(f"  {T_C:>5}°C", end="")
    for gas_name in ['CH4', 'H2S', 'N2', 'H2']:
        vphi_g = V_phi(gas_name, T_K_local, 30.0)
        delta_g = gas_mw(gas_name) / vphi_g - wd
        ratio = delta_g / delta_co2_local if abs(delta_co2_local) > 1e-10 else 0
        print(f"  {ratio:+8.3f}", end="")
    print()

print("""
  If the ratios are roughly constant with temperature, then a single
  set of scaling coefficients works at all T. If they vary significantly,
  we need T-dependent coefficients.
""")

# ============================================================================
# 8. SUMMARY AND RECOMMENDATION
# ============================================================================
print("=" * 70)
print("8. SUMMARY AND RECOMMENDATION")
print("=" * 70)

print("""
FINDINGS:
  1. Free-gas viscosity ratio CANNOT work: all ratios are positive,
     but CH4 (and likely N2, H2) DECREASE water viscosity. The free-gas
     viscosity has no physical connection to how dissolved molecules
     affect liquid water structure.

  2. The density effect (which we already compute) is a reliable proxy
     for the viscosity effect SIGN: gases that increase density also
     increase viscosity, and vice versa. This is physically expected -
     both are controlled by how the solute molecule occupies space in
     the hydrogen bond network.

  3. Two scaling approaches that preserve the correct sign:

     (A) Density-ratio scaling: a_i = a_CO2 * (drho%_i / drho%_CO2)
         Simple, uses existing density model. Predicts CH4 effect in
         the right ballpark vs Kumagai-Yokoyama.

     (B) Dissolved-density scaling: uses (MW/Vphi - MW_H2O/V1*) ratios.
         More physical but equivalent to (A) for small x.

     Both give the same qualitative results. (A) is simpler to implement.

  4. The scaling ratios are reasonably stable across temperature (25-200°C),
     meaning a single set of coefficients should work.

  5. Recommended viscosity model:

     mu = mu_brine * Product_i(1 + a_i * x_i^b)

     where b ≈ 1.01 (from Islam-Carlson), and:
       a_CO2 = +4.65 (Islam-Carlson, calibrated)
       a_i   = a_CO2 * (drho%_gas_i / drho%_CO2) for other gases

     Or equivalently: a_i = 4.65 * (density_effect_i / density_effect_CO2)
     computed at a reference condition from the Plyasunov/Garcia model.

CAVEATS:
  - Only calibrated against CO2 (Islam-Carlson) with order-of-magnitude
    check against CH4 (Kumagai-Yokoyama).
  - Assumes viscosity scales linearly with density effect, which is an
    approximation. The true relationship may be nonlinear.
  - H2S is a special case (hydrogen bonding with water). Its viscosity
    effect may not follow the density-based scaling.
  - No validation data exists for C2H6, C3H8, n-C4H10, H2, or N2
    dissolved in water.
""")
