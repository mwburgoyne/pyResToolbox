"""
Brine density correlations for NaCl solutions.

Provides:
    - rho_brine(T, P, S): brine density in kg/m³ (primary entry point)
    - rho_brine_spivey(T, P, S): Spivey et al. (modified) from McCain (2011)
    - rho_brine_bw(T, P, S): Batzle & Wang (1992) correlation (legacy)

Parameters:
    T: temperature in K
    P: pressure in MPa
    S: salinity as weight fraction of NaCl (e.g., 0.1 = 100,000 ppm = 10 wt%)

References:
    Spivey et al. (modified): McCain, Petroleum Reservoir Fluid Properties,
        Chapter 4, Eqs. 4.1-4.12. Multi-step thermodynamic approach with
        reference state at 70 MPa. Used by pyResToolbox.
    Batzle & Wang (1992): Geophysics 57(11), 1396-1408.
        Simpler correlation. Valid: 10-350°C, 5-100 MPa, 0-0.3 wt fraction.
"""

import numpy as np
from water_properties import rho_w

# =============================================================================
# Spivey coefficient tables (from pyResToolbox brine.py / McCain Ch.4)
# =============================================================================

# Rational function form (Eq 4.1): f(t) = (a1*(t/100)^2 + a2*(t/100) + a3) / (a4*(t/100)^2 + a5*(t/100) + 1)
# Input arrays are [0, a1, a2, a3, a4, a5] (index 0 unused)

_RHOW_T70 = [0, -0.127213, 0.645486, 1.03265, -0.070291, 0.639589]
_EWT = [0, 4.221, -3.478, 6.221, 0.5182, -0.4405]
_FWT = [0, -11.403, 29.932, 27.952, 0.20684, 0.3768]
_DM2T = [0, -0.00011149, 0.000175105, -0.00043766, 0, 0]
_DM32T = [0, -0.0008878, -0.0001388, -0.00296318, 0, 0.51103]
_DM1T = [0, 0.0021466, 0.012427, 0.042648, -0.081009, 0.525417]
_DM12T = [0, 0.0002356, -0.0003636, -0.0002278, 0, 0]
_EMT = [0, 0, 0, 0.1249, 0, 0]
_FM32T = [0, -0.617, -0.747, -0.4339, 0, 10.26]
_FM1T = [0, 0, 9.917, 5.1128, 0, 3.892]
_FM12T = [0, 0.0365, -0.0369, 0, 0, 0]


def _eq41(t_celsius, coeffs):
    """McCain Eq. 4.1 rational function. t in °C, coeffs = [0, a1..a5]."""
    t2 = t_celsius / 100.0
    return (coeffs[1] * t2**2 + coeffs[2] * t2 + coeffs[3]) / \
           (coeffs[4] * t2**2 + coeffs[5] * t2 + 1.0)


def rho_brine_spivey(T, P, S):
    """
    Brine density from Spivey et al. (modified) per McCain Ch. 4, Eqs 4.1-4.12.

    This is the same correlation used by pyResToolbox. It uses a thermodynamic
    reference state at 70 MPa and compressibility-based corrections to compute
    density at arbitrary (T, P).

    Parameters:
        T: temperature in K
        P: pressure in MPa
        S: salinity as weight fraction of NaCl (0 to ~0.3)

    Returns:
        brine density in kg/m³
    """
    T_C = float(T) - 273.15
    MPa = float(P)
    S_val = float(S)

    # Molality from weight fraction: m = 1000*S / (M_NaCl * (1-S))
    M_NaCl = 58.4428
    if S_val > 0:
        wt_pct = S_val  # weight fraction (0-1)
        m = 1000.0 * wt_pct / (M_NaCl * (1.0 - wt_pct))
    else:
        m = 0.0

    # Step 1: Density of pure water at T and 70 MPa reference pressure (g/cm³)
    rhow_t70 = _eq41(T_C, _RHOW_T70)

    # Step 2: Compressibility coefficients for pure water
    Ewt = _eq41(T_C, _EWT)
    Fwt = _eq41(T_C, _FWT)

    # Step 3: Pure water density at (T, P) via compressibility integral
    Iwt70 = (1.0 / Ewt) * np.log(abs(Ewt + Fwt))              # Eq 4.3
    Iwtp = (1.0 / Ewt) * np.log(abs(Ewt * (MPa / 70.0) + Fwt))  # Eq 4.4
    rhowtp = rhow_t70 * np.exp(Iwtp - Iwt70)                   # Eq 4.5

    if m == 0.0:
        return float(rhowtp) * 1000.0  # g/cm³ → kg/m³

    # Step 4: Salt correction coefficients
    Dm2t = _eq41(T_C, _DM2T)
    Dm32t = _eq41(T_C, _DM32T)
    Dm1t = _eq41(T_C, _DM1T)
    Dm12t = _eq41(T_C, _DM12T)

    # Step 5: Brine compressibility correction coefficients
    Emt = _eq41(T_C, _EMT)
    Fm32t = _eq41(T_C, _FM32T)
    Fm1t = _eq41(T_C, _FM1T)
    Fm12t = _eq41(T_C, _FM12T)

    # Brine density at T and 70 MPa (Eq 4.6)
    rhobt70 = rhow_t70 + Dm2t * m**2 + Dm32t * m**1.5 + Dm1t * m + Dm12t * m**0.5

    # Brine compressibility coefficients (Eqs 4.7-4.8)
    Ebtm = Ewt + Emt * m
    Fbtm = Fwt + Fm32t * m**1.5 + Fm1t * m + Fm12t * m**0.5

    # Brine density at (T, P) via compressibility integral (Eqs 4.10-4.12)
    Ibt70 = (1.0 / Ebtm) * np.log(abs(Ebtm + Fbtm))
    Ibtpm = (1.0 / Ebtm) * np.log(abs(Ebtm * (MPa / 70.0) + Fbtm))
    rhobtpm = rhobt70 * np.exp(Ibtpm - Ibt70)  # g/cm³

    return float(rhobtpm) * 1000.0  # g/cm³ → kg/m³


def rho_brine_bw(T, P, S):
    """
    Brine density from Batzle & Wang (1992) Eq. 27.

    Parameters:
        T: temperature in K (converted internally to °C)
        P: pressure in MPa
        S: salinity as weight fraction of NaCl (0 to ~0.3)

    Returns:
        brine density in kg/m³
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    P = np.atleast_1d(np.asarray(P, dtype=float))
    S = np.atleast_1d(np.asarray(S, dtype=float))

    T, P, S = np.broadcast_arrays(T, P, S)

    T_C = T - 273.15

    rho_water = np.zeros_like(T_C)
    for i in np.ndindex(T.shape):
        rho_water[i] = rho_w(float(T[i]), float(P[i]))

    rho_w_gcc = rho_water / 1000.0

    rho_b_gcc = rho_w_gcc + S * (0.668 + 0.44 * S + 1e-6 * (
        300.0 * P - 2400.0 * P * S
        + T_C * (80.0 + 3.0 * T_C - 3300.0 * S - 13.0 * P + 47.0 * P * S)
    ))

    rho_b = rho_b_gcc * 1000.0

    return float(rho_b) if rho_b.ndim == 0 else rho_b.squeeze()


def rho_brine(T, P, S):
    """
    Brine density — main entry point.

    Uses Spivey et al. (modified) correlation from McCain (2011).

    Parameters:
        T: temperature in K
        P: pressure in MPa
        S: salinity as weight fraction of NaCl

    Returns:
        density in kg/m³
    """
    return rho_brine_spivey(T, P, S)


def rho_water_spivey(T, P):
    """
    Pure water density from Spivey correlation (S=0).

    Useful for comparing against IAPWS-95.

    Parameters:
        T: temperature in K
        P: pressure in MPa

    Returns:
        density in kg/m³
    """
    return rho_brine_spivey(T, P, 0.0)


def salinity_from_ppm(ppm):
    """Convert salinity from ppm (mg/L ~ mg/kg) to weight fraction."""
    return ppm / 1e6


def salinity_from_molality(m_NaCl):
    """
    Convert NaCl molality (mol/kg solvent) to weight fraction.

    S = m * M_NaCl / (1000 + m * M_NaCl)
    where M_NaCl = 58.44 g/mol
    """
    M_NaCl = 58.44
    return m_NaCl * M_NaCl / (1000.0 + m_NaCl * M_NaCl)


if __name__ == "__main__":
    print("=== Brine Density Validation ===\n")

    T_test = 298.15  # 25°C
    P_test = 10.0    # 10 MPa

    rho_pure = rho_w(T_test, P_test)
    rho_pure_spivey = rho_water_spivey(T_test, P_test)
    print(f"Pure water at T={T_test}K, P={P_test}MPa:")
    print(f"  IAPWS-95: ρ = {rho_pure:.4f} kg/m³")
    print(f"  Spivey:   ρ = {rho_pure_spivey:.4f} kg/m³")
    print(f"  Diff:     {rho_pure_spivey - rho_pure:+.4f} kg/m³")

    print(f"\nBrine density (Spivey vs Batzle-Wang):")
    print(f"{'S(wt%)':>8} {'Spivey':>10} {'B&W':>10} {'Diff':>8} {'IAPWS_w':>10}")
    print("-" * 52)

    for S_pct in [0, 1, 3, 5, 10, 15, 20, 25]:
        S = S_pct / 100.0
        rho_sp = rho_brine_spivey(T_test, P_test, S)
        rho_bw = rho_brine_bw(T_test, P_test, S)
        print(f"{S_pct:8.1f} {rho_sp:10.4f} {rho_bw:10.4f} {rho_sp - rho_bw:+8.4f} {rho_pure:10.4f}")

    # Temperature sweep at 10 wt% NaCl
    S_10 = 0.10
    print(f"\n\nBrine density at S=10 wt% NaCl, P=30 MPa:")
    print(f"{'T(°C)':>6} {'Spivey':>10} {'B&W':>10} {'Diff':>8}")
    print("-" * 38)
    for T_C in [25, 50, 100, 150, 200, 250]:
        T_K = T_C + 273.15
        rho_sp = rho_brine_spivey(T_K, 30.0, S_10)
        rho_bw = rho_brine_bw(T_K, 30.0, S_10)
        print(f"{T_C:6} {rho_sp:10.4f} {rho_bw:10.4f} {rho_sp - rho_bw:+8.4f}")

    # Molality conversion test
    print(f"\n\nMolality to weight fraction:")
    for m in [0.5, 1.0, 2.0, 4.0]:
        S = salinity_from_molality(m)
        print(f"  {m:.1f} mol/kg → S = {S:.4f} ({S*100:.2f} wt%)")
