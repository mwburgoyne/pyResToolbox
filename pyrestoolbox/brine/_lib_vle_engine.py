"""
Unified VLE Engine for Multi-Gas Phase Equilibria
==================================================
Handles H2O + any combination of H2, CO2, N2, H2S, CH4, C2H6, C3H8, nC4H10,
iC4H10, iC5H12, nC5H12, nC6H14, nC7H16, nC8H18, nC10H22 with salting-out support.

SOLUTION SCHEME (confirmed by Curtis Whitson, Feb 2026):
Two independent full-mixture flashes using Rachford-Rice:
  Flash 1: All gas-water BIPs = kij_AQ → take AQUEOUS phase (gas solubilities)
  Flash 2: All gas-water BIPs = kij_NA → take NON-AQUEOUS phase (water content)
  True K-values: K_i = y_i(Flash 2) / x_i(Flash 1)

Gas-gas BIPs from literature are used in both flashes.

For DILUTE systems (H2, natural gas), the binary decomposition + Henry's law
gives equivalent results (tested: <0.4% for H2 in 90/10 H2/CH4 mix).
For CONCENTRATED systems (acid gas), the full multi-component flash is needed.

SALINITY: Three approaches supported:
  - Gamma-phi (RECOMMENDED): Freshwater EOS + Sechenov activity coefficient
    in K-value: K_i = γ_i × φ_i^L / φ_i^V, where γ_i = 10^(ks_i × m).
    Standard electrolyte thermodynamics. One multiplication in existing
    custom flash. For dilute systems, algebraically equivalent to explicit
    Sechenov (<0.3% difference for H2).
  - Explicit Sechenov: Freshwater flash + post-solve Eq 8 correction.
  - Embedded: kij(csw) correlations (Eqs 12-15) — original S&W for HC/CO2/N2.

The Sechenov coefficient ks source is gas-dependent:
  - Hydrocarbons & H2: S&W Equation 8 with Tb as correlating parameter
  - Other gases (H2S, etc.): May need gas-specific ks correlations
  The gamma-phi FRAMEWORK is universal; the ks SOURCE is gas-dependent.

RACHFORD-RICE SOLVER:
  Uses the method of Nielsen & Lia (2022), Fluid Phase Equilibria, which
  gracefully handles catastrophic numerical roundoff errors through a
  transformed variable approach.

References:
- Soreide & Whitson, Fluid Phase Equilibria 77 (1992) 217-240 + Errata
- Yan et al., Fluid Phase Equilibria 298 (2011) 180-189 (improved CO2)
- Nielsen & Lia, Fluid Phase Equilibria (2022) - Robust RR solver
- This work (H2 correlations): A=-14.59, B=2.184, C=0.365, kij_NA=0.468
- This work (CH4 rational): A=-2.1642, B=1.7325, C=0.2105 (MC-3 alpha)

Author: Mark Burgoyne, Markus H. Nielsen
Date: 2025-2026
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable
from dataclasses import dataclass


# =============================================================================
# Constants
# =============================================================================
R_GAS = 8.314462  # J/(mol·K)
OMEGA_A = 0.45724
OMEGA_B = 0.07780
MW_NACL = 58.44
MW_H2O = 18.015

# H2 critical temperature for BIP correlations (manuscript value, NIST)
# This is used ONLY for reduced temperature in kij_aq_h2 correlation
# NIST H2 critical temperature used consistently for both EOS and BIP correlations
BIP_TC_H2 = 33.145  # K (NIST)


# =============================================================================
# Component Properties Database
# =============================================================================
@dataclass
class ComponentProperties:
    """Critical properties and parameters for a component."""
    name: str
    Tc: float      # Critical temperature (K)
    Pc: float      # Critical pressure (Pa)
    omega: float   # Acentric factor
    Tb: float      # Normal boiling point (K) - for Sechenov
    MW: float      # Molecular weight (g/mol)


# Component database: All S&W 1992 Table 5 gases + H2 (this work)
# Tc, Pc, omega from S&W 1992; Tb from NIST for Sechenov calculations
COMPONENTS = {
    # Water
    'H2O': ComponentProperties('Water', 647.3, 22.12e6, 0.3434, 373.15, 18.015),

    # Hydrogen (this work) - Uses Tc=33.145 K (NIST) consistently for both EOS and BIP correlations
    'H2': ComponentProperties('Hydrogen', 33.145, 1.2964e6, -0.219, 20.3, 2.016),

    # Acid gases
    'CO2': ComponentProperties('Carbon Dioxide', 304.2, 7.38e6, 0.2273, 194.7, 44.01),
    'H2S': ComponentProperties('Hydrogen Sulfide', 373.2, 8.94e6, 0.1081, 212.8, 34.082),

    # Nitrogen
    'N2': ComponentProperties('Nitrogen', 126.1, 3.40e6, 0.0403, 77.36, 28.014),

    # Hydrocarbons (S&W 1992 Table 5)
    'CH4': ComponentProperties('Methane', 190.6, 4.60e6, 0.0108, 111.66, 16.043),
    'C2H6': ComponentProperties('Ethane', 305.4, 4.88e6, 0.0986, 184.6, 30.07),
    'C3H8': ComponentProperties('Propane', 369.8, 4.25e6, 0.1524, 231.1, 44.097),
    'iC4H10': ComponentProperties('i-Butane', 408.1, 3.65e6, 0.1770, 261.4, 58.123),
    'nC4H10': ComponentProperties('n-Butane', 425.2, 3.80e6, 0.1931, 272.7, 58.123),
    'iC5H12': ComponentProperties('i-Pentane', 460.4, 3.38e6, 0.2270, 301.0, 72.15),
    'nC5H12': ComponentProperties('n-Pentane', 469.6, 3.37e6, 0.2510, 309.2, 72.15),
    'nC6H14': ComponentProperties('n-Hexane', 507.4, 3.01e6, 0.2990, 341.9, 86.18),
    'nC7H16': ComponentProperties('n-Heptane', 540.3, 2.74e6, 0.3490, 371.6, 100.2),
    'nC8H18': ComponentProperties('n-Octane', 568.8, 2.49e6, 0.3980, 398.8, 114.2),
    'nC10H22': ComponentProperties('n-Decane', 617.7, 2.10e6, 0.4900, 447.3, 142.3),
}

# List of supported gas species (excludes water)
GAS_SPECIES = ['H2', 'CO2', 'N2', 'H2S', 'CH4', 'C2H6', 'C3H8',
               'iC4H10', 'nC4H10', 'iC5H12', 'nC5H12',
               'nC6H14', 'nC7H16', 'nC8H18', 'nC10H22']


# =============================================================================
# Alpha Functions
# =============================================================================
def alpha_water_soreide(Tr: float, salinity_molal: float = 0.0) -> float:
    """
    Soreide-Whitson modified alpha for water/brine (Equation 9).

    Args:
        Tr: Reduced temperature T/Tc_water
        salinity_molal: Salt concentration in mol/kg water

    Returns:
        Alpha parameter for water in PR EOS
    """
    Tr_safe = max(Tr, 0.01)
    sqrt_alpha = (1.0 + 0.4530 * (1.0 - Tr_safe * (1.0 - 0.0103 * salinity_molal**1.1))
                  + 0.0034 * (Tr_safe**(-3) - 1.0))
    return sqrt_alpha**2


def alpha_water_mc3(Tr: float) -> float:
    """
    Mathias-Copeman 3-parameter alpha for pure water.

    Refitted to IAPWS-95 vapor pressure data (Wagner & Pruss 2002) using
    S&W critical properties (Tc=647.3 K, Pc=221.2 bar).
    Pvap MARE = 0.006% over 20-200 deg C (vs 1.18% for S&W original).

    Args:
        Tr: Reduced temperature T/Tc_water (Tc=647.3 K)

    Returns:
        Alpha parameter for water in PR EOS
    """
    x = 1.0 - np.sqrt(max(Tr, 0.01))
    sqrt_alpha = 1.0 + 0.9110 * x - 0.2756 * x**2 + 0.2185 * x**3
    return sqrt_alpha**2


def alpha_standard_pr(Tr: float, omega: float) -> float:
    """
    Standard Peng-Robinson alpha function.

    Args:
        Tr: Reduced temperature T/Tc
        omega: Acentric factor

    Returns:
        Alpha parameter for PR EOS
    """
    m = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    return (1.0 + m * (1.0 - np.sqrt(max(Tr, 0.0))))**2


# =============================================================================
# kij_AQ Correlations (Aqueous Phase) - with errata corrections
# =============================================================================
def kij_aq_hydrocarbon(T_K: float, omega: float, Tc: float,
                       salinity_molal: float = 0.0) -> float:
    """
    Aqueous phase BIP for hydrocarbons (S&W 1992 Equations 11-12 with errata).

    General form for CH4, C2H6, C3H8, nC4H10, etc.

    Args:
        T_K: Temperature in Kelvin
        omega: Acentric factor of the hydrocarbon
        Tc: Critical temperature of the hydrocarbon (K)
        salinity_molal: Salt concentration in mol/kg water

    Returns:
        kij_AQ value
    """
    Tr = T_K / Tc
    cs = salinity_molal

    # Corrected constants from S&W errata
    A0 = 1.1120 - 1.7369 * omega**(-0.1)
    A1 = 1.1001 + 0.8360 * omega
    A2 = -0.15742 - 1.0988 * omega

    # Corrected salinity coefficients from errata
    alpha0 = 0.017407
    alpha1 = 0.033516
    alpha2 = 0.011478

    kij = (A0 * (1.0 + alpha0 * cs)
           + A1 * Tr * (1.0 + alpha1 * cs)
           + A2 * Tr**2 * (1.0 + alpha2 * cs))
    return kij


def kij_aq_ch4(T_K: float, salinity_molal: float = 0.0) -> float:
    """
    Aqueous phase BIP for CH4-water (this work).

    Rational form correlation:
        kij = (A + Tr) / (B + C*Tr)  where Tr = T/Tc, Tc = 190.60 K

    Parameters (L1-optimal fit to 109 freshwater pointwise kij values, T <= 200C,
    regressed with MC-3 alpha for water):
        A = -2.1642
        B = 1.7325
        C = 0.2105

    Valid range: 25-200C (298-473 K)

    SALINITY: Returns FRESHWATER kij only. The salinity_molal parameter
    is accepted for API consistency but is NOT used. Salinity is handled
    via gamma-phi (Sechenov) post-correction with S&W Eq 8.
    """
    Tr = T_K / 190.60
    A, B, C = -2.1642, 1.7325, 0.2105
    return (A + Tr) / (B + C * Tr)


def kij_aq_c2h6(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for C2H6-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['C2H6'].omega,
                              COMPONENTS['C2H6'].Tc, salinity_molal)


def kij_aq_c3h8(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for C3H8-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['C3H8'].omega,
                              COMPONENTS['C3H8'].Tc, salinity_molal)


def kij_aq_ic4h10(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for iC4H10-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['iC4H10'].omega,
                              COMPONENTS['iC4H10'].Tc, salinity_molal)


def kij_aq_nc4h10(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for nC4H10-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['nC4H10'].omega,
                              COMPONENTS['nC4H10'].Tc, salinity_molal)


def kij_aq_ic5h12(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for iC5H12-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['iC5H12'].omega,
                              COMPONENTS['iC5H12'].Tc, salinity_molal)


def kij_aq_nc5h12(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for nC5H12-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['nC5H12'].omega,
                              COMPONENTS['nC5H12'].Tc, salinity_molal)


def kij_aq_nc6h14(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for nC6H14-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['nC6H14'].omega,
                              COMPONENTS['nC6H14'].Tc, salinity_molal)


def kij_aq_nc7h16(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for nC7H16-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['nC7H16'].omega,
                              COMPONENTS['nC7H16'].Tc, salinity_molal)


def kij_aq_nc8h18(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for nC8H18-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['nC8H18'].omega,
                              COMPONENTS['nC8H18'].Tc, salinity_molal)


def kij_aq_nc10h22(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for nC10H22-water/brine."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['nC10H22'].omega,
                              COMPONENTS['nC10H22'].Tc, salinity_molal)


def kij_aq_co2(T_K: float, salinity_molal: float = 0.0) -> float:
    """
    Aqueous phase BIP for CO2-water/brine (S&W 1992 Equation 14).

    Note: Near-critical CO2 conditions (T < 40C) may show larger deviations.
    """
    cs = salinity_molal
    Tr = T_K / COMPONENTS['CO2'].Tc

    term1 = -0.31092 * (1.0 + 0.15587 * cs**0.7505)
    term2 = 0.23580 * (1.0 + 0.17837 * cs**0.979) * Tr
    term3 = -21.2566 * np.exp(-6.7222 * Tr - cs)

    return term1 + term2 + term3


def kij_aq_n2(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for N2-water/brine (S&W 1992 Equation 13, with errata)."""
    Tr = T_K / COMPONENTS['N2'].Tc
    cs = salinity_molal
    kij = (-1.70235 * (1.0 + 0.025587 * cs**0.75)
           + 0.44338 * (1.0 + 0.08126 * cs**0.75) * Tr)
    return kij


def kij_aq_h2s(T_K: float, salinity_molal: float = 0.0) -> float:
    """Aqueous phase BIP for H2S-water (S&W 1992 Equation 15)."""
    Tr = T_K / COMPONENTS['H2S'].Tc
    return -0.20441 + 0.23426 * Tr


def kij_aq_h2(T_K: float, salinity_molal: float = 0.0) -> float:
    """
    Aqueous phase BIP for H2-water (this work).

    Rational form correlation:
        kij = (A + Tr) / (B + C*Tr)

    Parameters (from manuscript, fitted to Wiebe+Chabab+Torin-Ollarves data):
        A = -14.59
        B = 2.184
        C = 0.365

    Valid range: 0-165C (273-438 K)

    NOTE: Uses BIP_TC_H2=33.145 K consistently for reduced temperature.

    SALINITY: This function returns the FRESHWATER kij only.
    The salinity_molal parameter is accepted for API consistency but
    is NOT used. For H2, salinity is handled via the gamma-phi
    approach: K_i = γ_i × φ_i^L / φ_i^V where γ_i = 10^(ks × m),
    with ks from S&W Equation 8 using H2 boiling point (20.3 K).
    """
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for correlation consistency
    A, B, C = -14.59, 2.184, 0.365
    return (A + Tr) / (B + C * Tr)


# Gases whose S&W original kij_AQ correlations have EMBEDDED salinity (csw terms).
# Used only in 'sw_original' framework mode.
_SW_GASES_WITH_EMBEDDED_SALINITY = {
    'C2H6', 'C3H8', 'iC4H10', 'nC4H10',
    'iC5H12', 'nC5H12', 'nC6H14', 'nC7H16', 'nC8H18', 'nC10H22',
    'CO2', 'N2',
}
# Backward-compatible alias
GASES_WITH_EMBEDDED_SALINITY = _SW_GASES_WITH_EMBEDDED_SALINITY


# --- Proposed freshwater kij_AQ (Paper 2, MC-3 alpha) ---
def kij_aq_co2_proposed(T_K: float, salinity_molal: float = 0.0) -> float:
    """CO2: cubic in T(K), n=614, MAE=0.0119 (MC-3 alpha)."""
    return -1.5127 + 9.4980e-3*T_K - 2.1680e-5*T_K**2 + 1.8500e-8*T_K**3

def kij_aq_h2s_proposed(T_K: float, salinity_molal: float = 0.0) -> float:
    """H2S: exp form, n=409, MAE=0.0138 (MC-3 alpha)."""
    return -70.8170/T_K + 1540.9516*np.exp(-4532.56/T_K) + 0.21517

def kij_aq_n2_proposed(T_K: float, salinity_molal: float = 0.0) -> float:
    """N2: linear in T(K), n=134, MAE=0.0163 (MC-3 alpha)."""
    return -1.6669 + 3.447873e-3 * T_K

def kij_aq_h2_proposed(T_K: float, salinity_molal: float = 0.0) -> float:
    """H2: rational, n=162, MAE=0.0400 (Paper 2 MC-3 coefficients)."""
    Tr = T_K / 33.145
    return (-14.6157 + Tr) / (3.5494 + 0.2230 * Tr)

def kij_aq_c2h6_proposed(T_K: float, salinity_molal: float = 0.0) -> float:
    """C2H6: rational, n=94, MAE=0.0095 (MC-3 alpha)."""
    Tr = T_K / 305.40
    return (-1.2685 + Tr) / (0.3647 + 1.2800 * Tr)

def kij_aq_c3h8_proposed(T_K: float, salinity_molal: float = 0.0) -> float:
    """C3H8: rational, n=59, MAE=0.0022 (MC-3 alpha)."""
    Tr = T_K / 369.80
    return (-1.1492 + Tr) / (0.6127 + 1.3198 * Tr)

# CH4 kij_aq_ch4 is already MC-3 rational form — used directly in proposed dict.
# For nC4H10+ in proposed mode: call kij_aq_hydrocarbon at cs=0 (freshwater).


def kij_aq_ch4_sw_original(T_K: float, salinity_molal: float = 0.0) -> float:
    """CH4: S&W Eqs 11-12 with embedded salinity (passthrough)."""
    return kij_aq_hydrocarbon(T_K, COMPONENTS['CH4'].omega,
                              COMPONENTS['CH4'].Tc, salinity_molal)


# --- Drop-in kij_AQ (Track 2, S&W original alpha, freshwater only) ---
def kij_aq_co2_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """CO2: cubic in T(K), n=611, MAE=0.0119 (S&W alpha)."""
    return -1.5893 + 9.8873e-3*T_K - 2.2188e-5*T_K**2 + 1.8499e-8*T_K**3

def kij_aq_h2s_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """H2S: exp form, n=405, MAE=0.0112 (S&W alpha)."""
    return -74.6914/T_K + 1348.9615*np.exp(-4504.96/T_K) + 0.22598

def kij_aq_n2_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """N2: linear in T(K), n=127, MAE=0.0114 (S&W alpha)."""
    return -1.6689 + 3.441589e-3 * T_K

def kij_aq_h2_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """H2: rational, n=154, MAE=0.0300 (S&W alpha)."""
    Tr = T_K / 33.145
    return (-14.9412 + Tr) / (2.2832 + 0.3893 * Tr)

def kij_aq_ch4_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """CH4: rational, n=115, MAE=0.0089 (S&W alpha)."""
    Tr = T_K / 190.60
    return (-2.1756 + Tr) / (1.0388 + 0.6436 * Tr)

def kij_aq_c2h6_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """C2H6: rational, n=94, MAE=0.0095 (S&W alpha)."""
    Tr = T_K / 305.40
    return (-1.2668 + Tr) / (0.1739 + 1.4165 * Tr)

def kij_aq_c3h8_dropin(T_K: float, salinity_molal: float = 0.0) -> float:
    """C3H8: rational, n=59, MAE=0.0022 (S&W alpha)."""
    Tr = T_K / 369.80
    return (-1.1496 + Tr) / (0.3501 + 1.5930 * Tr)


# --- Proposed freshwater wrappers for HCs without fitted correlations ---
def _kij_aq_hc_proposed(gas: str):
    """Return a freshwater-only wrapper around kij_aq_hydrocarbon for a given HC."""
    comp = COMPONENTS[gas]
    def _fn(T_K: float, salinity_molal: float = 0.0) -> float:
        return kij_aq_hydrocarbon(T_K, comp.omega, comp.Tc, 0.0)
    _fn.__doc__ = f"{gas}: S&W Eqs 11-12 at cs=0 (no fitted proposed form)."
    return _fn


# =============================================================================
# Dual Dispatch Dictionaries
# =============================================================================
# Proposed: MC-3 alpha + freshwater-only kij. Salinity via Sechenov/gamma-phi.
KIJ_AQ_PROPOSED: Dict[str, Callable] = {
    'H2': kij_aq_h2_proposed,
    'CO2': kij_aq_co2_proposed,
    'N2': kij_aq_n2_proposed,
    'H2S': kij_aq_h2s_proposed,
    'CH4': kij_aq_ch4,              # Already MC-3 rational form
    'C2H6': kij_aq_c2h6_proposed,
    'C3H8': kij_aq_c3h8_proposed,
    'iC4H10': _kij_aq_hc_proposed('iC4H10'),
    'nC4H10': _kij_aq_hc_proposed('nC4H10'),
    'iC5H12': _kij_aq_hc_proposed('iC5H12'),
    'nC5H12': _kij_aq_hc_proposed('nC5H12'),
    'nC6H14': _kij_aq_hc_proposed('nC6H14'),
    'nC7H16': _kij_aq_hc_proposed('nC7H16'),
    'nC8H18': _kij_aq_hc_proposed('nC8H18'),
    'nC10H22': _kij_aq_hc_proposed('nC10H22'),
}

# S&W original: S&W alpha + embedded salinity in kij for HC/CO2/N2.
KIJ_AQ_SW_ORIGINAL: Dict[str, Callable] = {
    'H2': kij_aq_h2,                # Paper 1 rational (S&W had no H2)
    'CO2': kij_aq_co2,              # S&W Eq 14 (embedded salinity)
    'N2': kij_aq_n2,                # S&W Eq 13 (embedded salinity)
    'H2S': kij_aq_h2s,              # S&W Eq 15 (no salinity)
    'CH4': kij_aq_ch4_sw_original,  # S&W Eqs 11-12 (embedded salinity)
    'C2H6': kij_aq_c2h6,            # S&W Eqs 11-12 (embedded salinity)
    'C3H8': kij_aq_c3h8,            # S&W Eqs 11-12 (embedded salinity)
    'iC4H10': kij_aq_ic4h10,        # S&W Eqs 11-12 (embedded salinity)
    'nC4H10': kij_aq_nc4h10,        # S&W Eqs 11-12 (embedded salinity)
    'iC5H12': kij_aq_ic5h12,        # S&W Eqs 11-12 (embedded salinity)
    'nC5H12': kij_aq_nc5h12,        # S&W Eqs 11-12 (embedded salinity)
    'nC6H14': kij_aq_nc6h14,        # S&W Eqs 11-12 (embedded salinity)
    'nC7H16': kij_aq_nc7h16,        # S&W Eqs 11-12 (embedded salinity)
    'nC8H18': kij_aq_nc8h18,        # S&W Eqs 11-12 (embedded salinity)
    'nC10H22': kij_aq_nc10h22,      # S&W Eqs 11-12 (embedded salinity)
}

# Drop-in: S&W alpha + freshwater-only kij refitted for S&W alpha + embedded delta.
KIJ_AQ_DROPIN: Dict[str, Callable] = {
    'H2': kij_aq_h2_dropin,
    'CO2': kij_aq_co2_dropin,
    'N2': kij_aq_n2_dropin,
    'H2S': kij_aq_h2s_dropin,
    'CH4': kij_aq_ch4_dropin,
    'C2H6': kij_aq_c2h6_dropin,
    'C3H8': kij_aq_c3h8_dropin,
    'iC4H10': _kij_aq_hc_proposed('iC4H10'),
    'nC4H10': _kij_aq_hc_proposed('nC4H10'),
    'iC5H12': _kij_aq_hc_proposed('iC5H12'),
    'nC5H12': _kij_aq_hc_proposed('nC5H12'),
    'nC6H14': _kij_aq_hc_proposed('nC6H14'),
    'nC7H16': _kij_aq_hc_proposed('nC7H16'),
    'nC8H18': _kij_aq_hc_proposed('nC8H18'),
    'nC10H22': _kij_aq_hc_proposed('nC10H22'),
}

# Default dispatch = proposed (backward compatible alias)
KIJ_AQ_FUNCTIONS: Dict[str, Callable] = KIJ_AQ_PROPOSED


# =============================================================================
# Embedded Salinity BIP Parameters (Paper 2, all 8 gases)
# =============================================================================
# Form: kij(T,m) = kij_fw(T) + (a0 + a1*Tr + a2*Tr^2)*m  [+ (b0 + b1*Tr)*m^2 for CO2]
# Tr = T/Tc for each gas. Fitted with MC-3 alpha + proposed kij_fw.
EMBEDDED_SALINITY_PARAMS = {
    'CO2':    {'Tc': 304.20, 'a0': 0.0305, 'a1': -0.0525, 'a2': 0.0315,
               'b0': 0.0024, 'b1': -0.0027},  # Quadratic-in-m (5 params)
    'H2S':    {'Tc': 373.20, 'a0': 0.0454, 'a1': -0.0884, 'a2': 0.0484},
    'CH4':    {'Tc': 190.60, 'a0': 0.1383, 'a1': -0.1394, 'a2': 0.0442},
    'N2':     {'Tc': 126.10, 'a0': 0.2257, 'a1': -0.1548, 'a2': 0.0321},
    'H2':     {'Tc': 33.145, 'a0': 0.3833, 'a1': -0.0660, 'a2': 0.0033},
    'C2H6':   {'Tc': 305.40, 'a0': 0.0791, 'a1': -0.1263, 'a2': 0.0677},
    'C3H8':   {'Tc': 369.80, 'a0': 0.0606, 'a1': -0.1165, 'a2': 0.0772},
    'nC4H10': {'Tc': 425.20, 'a0': 0.0488, 'a1': -0.1072, 'a2': 0.0836},
}

# Drop-in embedded salinity BIP parameters (Track 2, S&W alpha)
# Placeholder — will be populated after running 07-Fit_Embedded_Salinity_All_Gases.py
# with framework='dropin'. These are fitted with brine VLE using alpha_water_soreide(Tr, m).
EMBEDDED_SALINITY_PARAMS_DROPIN: Dict[str, Dict] = {
    'CO2':  {'Tc': 304.20, 'a0': 0.0409, 'a1': -0.0807, 'a2': 0.0526,
             'b0': 0.0079, 'b1': -0.0085},  # Quadratic-in-m (5 params)
    'H2S':  {'Tc': 373.20, 'a0': 0.0341, 'a1': -0.0655, 'a2': 0.0376},
    'CH4':  {'Tc': 190.60, 'a0': 0.1304, 'a1': -0.1295, 'a2': 0.0394},
    'N2':   {'Tc': 126.10, 'a0': 0.2173, 'a1': -0.1468, 'a2': 0.0302},
    'H2':   {'Tc': 33.145, 'a0': 0.3658, 'a1': -0.0625, 'a2': 0.0030},
    'C2H6':   {'Tc': 305.40, 'a0': 0.0812, 'a1': -0.1286, 'a2': 0.0646},
    'C3H8':   {'Tc': 369.80, 'a0': 0.0606, 'a1': -0.1165, 'a2': 0.0772},  # From proposed (not refitted)
    'nC4H10': {'Tc': 425.20, 'a0': 0.0488, 'a1': -0.1072, 'a2': 0.0836},  # From proposed (not refitted)
}

# Paper 1 H2 embedded BIP coefficients (for sw_original mode)
_PAPER1_H2_EMBEDDED = {'Tc': 33.145, 'a0': 0.3833, 'a1': -0.06595, 'a2': 0.003321}


def calc_embedded_delta_kij(gas: str, T_K: float, salinity_molal: float,
                             params: Optional[Dict] = None) -> float:
    """
    Calculate the salinity correction delta_kij for embedded BIP approach.

    kij(T,m) = kij_fw(T) + delta_kij(T,m)
    delta_kij = (a0 + a1*Tr + a2*Tr^2)*m  [+ (b0 + b1*Tr)*m^2 for CO2]

    Args:
        gas: Gas name
        T_K: Temperature in Kelvin
        salinity_molal: NaCl molality
        params: Override parameter dict (default: use EMBEDDED_SALINITY_PARAMS)

    Returns:
        delta_kij value
    """
    if salinity_molal <= 0:
        return 0.0
    if params is None:
        if gas not in EMBEDDED_SALINITY_PARAMS:
            return 0.0
        params = EMBEDDED_SALINITY_PARAMS[gas]

    Tr = T_K / params['Tc']
    m = salinity_molal
    delta = (params['a0'] + params['a1'] * Tr + params['a2'] * Tr**2) * m
    if 'b0' in params:
        delta += (params['b0'] + params['b1'] * Tr) * m**2
    return delta


def get_kij_aq(gas: str, T_K: float, salinity_molal: float = 0.0,
               framework: str = 'proposed') -> float:
    """
    Get kij_AQ for any supported gas.

    Args:
        gas: Gas name (e.g., 'H2', 'CO2', 'CH4')
        T_K: Temperature in Kelvin
        salinity_molal: Salt concentration in mol/kg water
        framework: 'proposed' (MC-3 alpha, freshwater kij + Sechenov),
                   'sw_original' (S&W alpha with embedded salinity in kij), or
                   'dropin' (S&W alpha + new freshwater kij + embedded delta)

    Returns:
        kij_AQ value
    """
    if framework == 'proposed':
        dispatch = KIJ_AQ_PROPOSED
        if gas not in dispatch:
            raise ValueError(f"Unknown gas: {gas}. Supported: {list(dispatch.keys())}")
        # Proposed mode: always freshwater kij (salinity handled via Sechenov)
        return dispatch[gas](T_K, 0.0)
    elif framework == 'dropin':
        dispatch = KIJ_AQ_DROPIN
        if gas not in dispatch:
            raise ValueError(f"Unknown gas: {gas}. Supported: {list(dispatch.keys())}")
        # Drop-in mode: freshwater kij only (salinity via embedded delta)
        return dispatch[gas](T_K, 0.0)
    elif framework == 'sw_original':
        dispatch = KIJ_AQ_SW_ORIGINAL
        if gas not in dispatch:
            raise ValueError(f"Unknown gas: {gas}. Supported: {list(dispatch.keys())}")
        # S&W mode: pass salinity for embedded gases
        return dispatch[gas](T_K, salinity_molal)
    else:
        raise ValueError(f"Unknown framework: {framework}. "
                         "Use 'proposed', 'sw_original', or 'dropin'")


# =============================================================================
# kij_NA Values (Non-Aqueous Phase) - S&W 1992 Table 5 + this work
# =============================================================================
# Constant kij_NA values from S&W 1992 Table 5
KIJ_NA: Dict[str, Optional[float]] = {
    'H2': 0.468,       # This work
    'CO2': 0.1896,     # S&W 1992 (or 0.18756 from Yan 2011)
    'N2': 0.4778,      # S&W 1992
    'H2S': 0.1610,     # This work (constant; replaces S&W Eq 17)
    'CH4': 0.4850,     # S&W 1992
    'C2H6': 0.4920,    # S&W 1992
    'C3H8': 0.5070,    # S&W 1992
    'iC4H10': 0.5080,  # S&W 1992
    'nC4H10': 0.5080,  # S&W 1992
    'iC5H12': 0.5090,  # S&W 1992
    'nC5H12': 0.5090,  # S&W 1992
    'nC6H14': 0.5100,  # S&W 1992
    'nC7H16': 0.5100,  # S&W 1992
    'nC8H18': 0.5100,  # S&W 1992
    'nC10H22': 0.5100, # S&W 1992
}


def kij_na_h2s(T_K: float) -> float:
    """Non-aqueous phase BIP for H2S-water.

    Returns constant 0.161 (this work, optimized against y_H2O MARE).
    Replaces S&W 1992 Equation 17 (0.19031 - 0.05965*Tr).
    """
    return 0.1610


def kij_na_h2s_sw_eq17(T_K: float) -> float:
    """Original S&W 1992 Equation 17 for H2S kij_NA (retained for comparison)."""
    Tr = T_K / COMPONENTS['H2S'].Tc
    return 0.19031 - 0.05965 * Tr


def get_kij_na(gas: str, T_K: float) -> float:
    """
    Get kij_NA for any supported gas.

    Args:
        gas: Gas name
        T_K: Temperature in Kelvin

    Returns:
        kij_NA value
    """
    if gas not in KIJ_NA:
        raise ValueError(f"Unknown gas: {gas}. Supported: {list(KIJ_NA.keys())}")
    val = KIJ_NA.get(gas)
    if val is None:
        raise ValueError(f"kij_NA not defined for {gas}")
    return val


# =============================================================================
# Sechenov Coefficient (S&W 1992 Equation 8)
# =============================================================================
def sw_equation_8_ks(T_C: float, Tb_K: float) -> float:
    """
    Soreide-Whitson (1992) Equation 8 - Sechenov coefficient correlation.

    ks = 0.13163 + 4.45e-4*Tb - 7.692e-4*T + 2.6614e-6*T^2 - 2.612e-9*T^3

    where T is in F and Tb is normal boiling point in K.

    NOTE: Returns ks on log10 basis (per S&W Eq. 2).

    Args:
        T_C: Temperature in Celsius
        Tb_K: Normal boiling point of gas in Kelvin

    Returns:
        Sechenov coefficient ks (kg/mol, log10 basis)

    Example:
        For H2: Tb = 20.3 K
        For CH4: Tb = 111.66 K
        For N2: Tb = 77.36 K
    """
    T_F = T_C * 9/5 + 32
    ks = (0.13163 + 4.45e-4 * Tb_K - 7.692e-4 * T_F
          + 2.6614e-6 * T_F**2 - 2.612e-9 * T_F**3)
    return ks


def get_sechenov_ks(gas: str, T_K: float, salinity_molal: float = 1.0,
                    P_bar: float = 100.0,
                    framework: str = 'proposed') -> float:
    """
    Get the recommended Sechenov coefficient ks for a given gas.

    In 'proposed' framework, routes to best-available ks model per gas:
      - CO2: Dubessy et al. 2005 extended Sechenov + constant offset
      - H2S: Akinfiev et al. 2016 Pitzer model (salting_library)
      - All others (HCs, N2, H2): S&W Equation 8 (Tb-based)

    In 'sw_original' framework, all gases use S&W Equation 8.

    All return ks on log10 basis with molality scale.

    Args:
        gas: Gas name (e.g., 'H2', 'CO2', 'CH4')
        T_K: Temperature in Kelvin
        salinity_molal: NaCl molality (needed for Pitzer models). Default 1.0.
        P_bar: Pressure in bar. Default 100.
        framework: 'proposed', 'sw_original', or 'dropin'

    Returns:
        ks: Sechenov coefficient (log10 basis, kg/mol)
    """
    if framework == 'proposed':
        if gas == 'CO2':
            from pyrestoolbox.brine._lib_salting_library import ks_dubessy_co2
            # Dubessy 2005 + MARE-optimal shift (-0.011)
            return float(ks_dubessy_co2(T_K, m_NaCl=salinity_molal)) - 0.011
        elif gas == 'H2S':
            from pyrestoolbox.brine._lib_salting_library import ks_akinfiev_h2s
            # Akinfiev 2016 + MARE-optimal shift (+0.019)
            return float(ks_akinfiev_h2s(T_K, m_NaCl=salinity_molal)) + 0.019

    # All others, or sw_original for all gases: S&W Equation 8 with Tb
    if gas not in COMPONENTS:
        raise ValueError(f"Unknown gas: {gas}")
    T_C = T_K - 273.15
    ks = sw_equation_8_ks(T_C, COMPONENTS[gas].Tb)
    # N2: empirical +0.02 offset (proposed mode only)
    if gas == 'N2' and framework == 'proposed':
        ks += 0.02
    return ks


def salting_out_factor(T_C: float, Tb_K: float, salinity_molal: float) -> float:
    """
    Calculate salting-out factor from Sechenov coefficient.

    S&W Eq. 2: log10(x_fresh/x_brine) = ks * m
    Therefore: x_brine/x_fresh = 10^(-ks * m)

    Args:
        T_C: Temperature in Celsius
        Tb_K: Normal boiling point of gas in Kelvin
        salinity_molal: Salt concentration in mol/kg water

    Returns:
        Salting-out factor: x_brine/x_fresh (< 1 means reduced solubility)
    """
    if salinity_molal <= 0:
        return 1.0
    ks = sw_equation_8_ks(T_C, Tb_K)
    return 10**(-ks * salinity_molal)


# =============================================================================
# EOS Core Functions
# =============================================================================
@dataclass
class CubicSolution:
    """Extended cubic EOS solution with Gibbs energy analysis."""
    roots: List[float]
    Z_liquid: float
    Z_vapor: float
    delta_G: Optional[float]  # G_vapor - G_liquid (dimensionless, if two roots)
    two_roots: bool
    preferred_phase: str      # 'liquid' or 'vapor' (min Gibbs energy)


def _halley_cubic(c2: float, c1: float, c0: float, B: float) -> List[float]:
    """
    Solve depressed cubic Z^3 + c2*Z^2 + c1*Z + c0 = 0 using Halley iteration.

    Michelsen-style: start from inflection point, find largest root via Halley,
    then synthetic division + quadratic for remaining roots.

    Returns sorted list of valid roots (Z > B).
    """
    MAX_ITER = 50
    TOL = 1e-10

    # Inflection point: F''=0 at Z_inf = -c2/3
    Z_inf = -c2 / 3.0

    # Evaluate F at inflection
    F_inf = Z_inf**3 + c2 * Z_inf**2 + c1 * Z_inf + c0

    # Choose starting point for largest root
    if F_inf > 0:
        # Function positive at inflection — largest root is above inflection
        # Start above inflection, try B+1 or further right
        Z = max(B + 1.0, Z_inf + 1.0)
    else:
        # F_inf <= 0: largest root might be near or above local max
        # Discriminant of F': 2*c2^2 - 6*c1 = 2*(c2^2 - 3*c1)
        disc_Fp = c2**2 - 3.0 * c1
        if disc_Fp > 0:
            # Local max at Z = (-c2 + sqrt(disc)) / 3
            Z_local_max = (-c2 + np.sqrt(disc_Fp)) / 3.0
            F_max = Z_local_max**3 + c2 * Z_local_max**2 + c1 * Z_local_max + c0
            if F_max > 0:
                # Three real roots — start above local max for largest
                Z = Z_local_max + 0.5
            else:
                # Only one real root — start from B+1
                Z = max(B + 1.0, Z_inf + 1.0)
        else:
            # Monotonic cubic — one real root, start right of inflection
            Z = max(B + 1.0, Z_inf + 1.0)

    # Halley iteration for largest root
    converged = False
    for _ in range(MAX_ITER):
        F = Z**3 + c2 * Z**2 + c1 * Z + c0
        Fp = 3.0 * Z**2 + 2.0 * c2 * Z + c1
        Fpp = 6.0 * Z + 2.0 * c2

        if abs(Fp) < 1e-30:
            break

        DZ = F / Fp
        # Halley correction
        denom = 1.0 - 0.5 * DZ * Fpp / Fp
        if abs(denom) > 1e-15:
            DZ = DZ / denom

        Z -= DZ

        if abs(DZ) < TOL:
            converged = True
            break

    if not converged:
        return []  # Signal fallback needed

    Z1 = Z

    # Synthetic division: Z^3 + c2*Z^2 + c1*Z + c0 = (Z - Z1)(Z^2 + q1*Z + q0)
    q1 = c2 + Z1
    q0 = c1 + Z1 * q1

    # Solve quadratic Z^2 + q1*Z + q0 = 0
    disc = q1**2 - 4.0 * q0
    roots = [Z1]

    if disc >= 0:
        sqrt_disc = np.sqrt(disc)
        Z2 = (-q1 - sqrt_disc) / 2.0
        Z3 = (-q1 + sqrt_disc) / 2.0

        # Refine each quadratic root with one Halley step
        for Zk in [Z2, Z3]:
            F = Zk**3 + c2 * Zk**2 + c1 * Zk + c0
            Fp = 3.0 * Zk**2 + 2.0 * c2 * Zk + c1
            Fpp = 6.0 * Zk + 2.0 * c2
            if abs(Fp) > 1e-30:
                DZ = F / Fp
                denom = 1.0 - 0.5 * DZ * Fpp / Fp
                if abs(denom) > 1e-15:
                    DZ = DZ / denom
                Zk -= DZ
            roots.append(Zk)

    # Filter valid roots: Z > B
    valid = [r for r in roots if r > B + 1e-10]
    return sorted(valid)


def _gibbs_delta(A: float, B: float, Zl: float, Zv: float) -> float:
    """
    Gibbs energy difference G_vapor - G_liquid (dimensionless) for PR EOS.

    Michelsen formula:
    ΔG = ln((Zv-B)/(Zl-B)) - A/(2√2·B)·ln((Zv+δ₂B)(Zl+δ₁B)/((Zv+δ₁B)(Zl+δ₂B))) - (Zv-Zl)
    where δ₁ = 1+√2, δ₂ = 1-√2.

    If ΔG < 0, vapor is preferred; if ΔG > 0, liquid is preferred.
    """
    SQRT2 = 1.4142135623730951
    d1 = 1.0 + SQRT2  # δ₁
    d2 = 1.0 - SQRT2  # δ₂

    if Zv - B <= 0 or Zl - B <= 0 or B < 1e-15:
        return 0.0

    term1 = np.log((Zv - B) / (Zl - B))

    num = (Zv + d2 * B) * (Zl + d1 * B)
    den = (Zv + d1 * B) * (Zl + d2 * B)
    if den <= 0 or num <= 0:
        return 0.0

    term2 = -(A / (2.0 * SQRT2 * B)) * np.log(num / den)
    term3 = -(Zv - Zl)

    return term1 + term2 + term3


def solve_cubic_eos(A: float, B: float, return_info: bool = False):
    """
    Solve PR cubic EOS for compressibility factor Z.

    Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0

    Uses Halley-accelerated iteration (Michelsen-style) with np.roots fallback.

    Args:
        A: EOS parameter a*P/(R*T)^2
        B: EOS parameter b*P/(R*T)
        return_info: If True, return CubicSolution with Gibbs energy analysis

    Returns:
        List of valid Z roots (sorted ascending), or CubicSolution if return_info=True
    """
    c2 = -(1.0 - B)
    c1 = A - 3.0 * B**2 - 2.0 * B
    c0 = -(A * B - B**2 - B**3)

    # Try Halley solver first
    valid = _halley_cubic(c2, c1, c0, B)

    # Fallback to np.roots if Halley failed
    if not valid:
        roots = np.roots([1.0, c2, c1, c0])
        valid = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > B + 1e-10]
        valid = sorted(valid) if valid else [max(B + 0.01, 0.1)]

    if not return_info:
        return valid

    # Build CubicSolution with Gibbs energy analysis
    Zl = valid[0]
    Zv = valid[-1]
    two_roots = len(valid) >= 2 and abs(Zv - Zl) > 1e-10

    if two_roots:
        dG = _gibbs_delta(A, B, Zl, Zv)
        preferred = 'vapor' if dG < 0 else 'liquid'
    else:
        dG = None
        # Single root — classify by proximity to B (liquid-like) vs 1 (vapor-like)
        preferred = 'liquid' if Zl < 0.5 else 'vapor'

    return CubicSolution(
        roots=valid,
        Z_liquid=Zl,
        Z_vapor=Zv,
        delta_G=dG,
        two_roots=two_roots,
        preferred_phase=preferred,
    )


def calc_fugacity_coeff(Z: float, A: float, B: float,
                        Bi_over_B: float, sum_xAij_over_A: float) -> float:
    """
    Calculate fugacity coefficient for component i in mixture.

    Args:
        Z: Compressibility factor
        A: Mixture A parameter
        B: Mixture B parameter
        Bi_over_B: Bi/B ratio for component i
        sum_xAij_over_A: (sum_j x_j * Aij) / A for component i

    Returns:
        Fugacity coefficient phi_i
    """
    sqrt2 = np.sqrt(2.0)

    if B < 1e-15 or A < 1e-15 or Z <= B:
        return 1.0

    term1 = Bi_over_B * (Z - 1.0)
    term2 = -np.log(max(Z - B, 1e-15))

    log_arg = (Z + (1.0 + sqrt2) * B) / (Z + (1.0 - sqrt2) * B)
    if log_arg > 0:
        term3 = (A / (2.0 * sqrt2 * B)) * (Bi_over_B - 2.0 * sum_xAij_over_A) * np.log(log_arg)
    else:
        term3 = 0.0

    ln_phi = term1 + term2 + term3
    return np.exp(min(ln_phi, 50))


# =============================================================================
# Rachford-Rice Solver — Nielsen & Lia (2022), Fluid Phase Equilibria
# =============================================================================
def rr_solver(
    zi: np.ndarray, ki: np.ndarray,
    tol: float = 1e-15, max_iter: int = 100
) -> Tuple[int, np.ndarray, np.ndarray, float, float]:
    """
    Solve the Rachford-Rice equation using the method of Nielsen & Lia (2022),
    which gracefully handles catastrophic numerical roundoff errors through
    a transformed variable approach.

    Reference:
        M. Nielsen & H. Lia, "Generalized Rachford-Rice Algorithm",
        Fluid Phase Equilibria (2022)

    Args:
        zi: Molar composition (will be normalized)
        ki: K-values for each component
        tol: Solution tolerance (default 1e-15)
        max_iter: Maximum iterations (default 100)

    Returns:
        N_it: Number of iterations required
        yi: Vapor mole fraction compositions
        xi: Liquid mole fraction compositions
        V: Vapor molar fraction
        L: Liquid molar fraction
    """
    zi = zi / np.sum(zi)  # Normalize feed compositions

    # Perturb K-values exactly equal to 1.0 to avoid singularity in ci = 1/(1-K)
    # Components with K=1 contribute nothing to RR equation; perturbation is safe
    ki = np.where(np.abs(ki - 1.0) < 1e-12, 1.0 + 1e-12, ki)

    def rr(V: float) -> float:
        return np.dot(zi, (ki - 1) / (1 + V * (ki - 1)))

    # Assess if solution is nearest vapor or liquid
    near_vapor = rr(0.5) > 0

    ki_hat = 1.0 / ki if near_vapor else ki.copy()
    ci = 1.0 / (1.0 - ki_hat)                       # Eq 10

    # Transformed variable bounds
    phi_max = min(1.0 / (1.0 - np.min(ki_hat)), 0.5)  # Eq 11a
    phi_min = 1.0 / (1.0 - np.max(ki_hat))             # Eq 11b
    b_min = 1.0 / (phi_max - phi_min)                   # Eq 15
    b_max = np.inf

    b = 1.0 / (0.25 - phi_min)

    def h(b: float) -> float:                          # Eq 12b
        return np.sum(zi * b / (1.0 + b * (phi_min - ci)))

    def dh(b: float) -> float:                         # Eq 16b
        return np.sum(zi / (1.0 + b * (phi_min - ci))**2)

    N_it = 0
    h_b = np.inf

    while abs(h_b) > tol:
        N_it += 1
        h_b = h(b)
        dh_b = dh(b)

        if h_b > 0:
            b_max = b
        else:
            b_min = b

        b = b - h_b / dh_b

        if b < b_min or b > b_max:
            b = (b_min + b_max) / 2.0

        if N_it > max_iter:
            break

    # Recover compositions from transformed variables
    ui = -zi * ci * b / (1.0 + b * (phi_min - ci))    # Eq 27b
    phi = (1.0 + b * phi_min) / b                      # Rearranged Eq 14b

    if near_vapor:
        L = phi
        V = 1.0 - L
        yi = ui
        xi = ki_hat * ui                                # Eq 28
    else:
        V = phi
        L = 1.0 - V
        xi = ui
        yi = ki_hat * ui                                # Eq 28

    return N_it, yi, xi, V, L


def solve_rachford_rice(z: np.ndarray, K: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Convenience wrapper around rr_solver for flash calculations.

    Handles single-phase detection before calling the robust solver.

    Args:
        z: Feed composition (will be normalized)
        K: K-values

    Returns:
        V: Vapor fraction
        x: Liquid mole fractions
        y: Vapor mole fractions
    """
    z = np.asarray(z, dtype=float)
    K = np.asarray(K, dtype=float)
    z = z / np.sum(z)

    Km1 = K - 1.0

    # Single-phase checks
    if np.sum(z * Km1) <= 0:
        # All liquid — V = 0
        return 0.0, z.copy(), (K * z) / np.sum(K * z)
    if np.sum(z * Km1 / K) >= 0:
        # All vapor — V = 1
        return 1.0, (z / K) / np.sum(z / K), z.copy()

    # Two-phase: use robust solver
    N_it, yi, xi, V, L = rr_solver(z, K)
    return V, xi, yi


# =============================================================================
# Unit Conversions
# =============================================================================
def fahrenheit_to_kelvin(T_F: float) -> float:
    """Convert Fahrenheit to Kelvin."""
    return (T_F - 32) * 5/9 + 273.15


def celsius_to_kelvin(T_C: float) -> float:
    """Convert Celsius to Kelvin."""
    return T_C + 273.15


def kelvin_to_celsius(T_K: float) -> float:
    """Convert Kelvin to Celsius."""
    return T_K - 273.15


def psia_to_pascal(P_psia: float) -> float:
    """Convert psia to Pascal."""
    return P_psia * 6894.757


def bar_to_pascal(P_bar: float) -> float:
    """Convert bar to Pascal."""
    return P_bar * 1e5


def pascal_to_bar(P_Pa: float) -> float:
    """Convert Pascal to bar."""
    return P_Pa / 1e5


def wt_pct_to_molality(wt_pct_NaCl: float) -> float:
    """
    Convert weight percent NaCl to molality.

    Args:
        wt_pct_NaCl: Salt concentration in weight percent

    Returns:
        Molality in mol NaCl / kg water
    """
    if wt_pct_NaCl <= 0:
        return 0.0
    if wt_pct_NaCl >= 100:
        raise ValueError("Salt concentration cannot be 100% or more")
    mass_salt = wt_pct_NaCl
    mass_water = 100 - wt_pct_NaCl
    moles_salt = mass_salt / MW_NACL
    kg_water = mass_water / 1000
    return moles_salt / kg_water


def molality_to_wt_pct(molality: float) -> float:
    """
    Convert molality to weight percent NaCl.

    Args:
        molality: mol NaCl / kg water

    Returns:
        Weight percent NaCl
    """
    if molality <= 0:
        return 0.0
    mass_salt = molality * MW_NACL  # grams per kg water
    mass_total = 1000 + mass_salt   # total grams
    return mass_salt / mass_total * 100


def y_h2o_to_stb_mmscf(y_H2O: float) -> float:
    """
    Convert water mole fraction in wet gas to stb/mmscf of dry gas.

    Derivation:
      - Standard conditions: 60F, 14.7 psia
      - Molar volume at SC: 379.38 scf/lbmol
      - Water MW: 18.015 lb/lbmol
      - Water density at 60F: 350.2 lb/bbl

      stb/mmscf = [y_H2O/(1-y_H2O)] x [1e6/379.38] x [18.015/350.2]
                = [y_H2O/(1-y_H2O)] x 135.61

    Args:
        y_H2O: Water mole fraction in wet gas (mol/mol)

    Returns:
        Water content in stb per mmscf of dry gas
    """
    if y_H2O <= 0:
        return 0.0
    if y_H2O >= 1:
        raise ValueError("y_H2O must be less than 1")

    CONV_FACTOR = 135.61
    return y_H2O / (1.0 - y_H2O) * CONV_FACTOR


def y_h2o_to_lb_mmscf(y_H2O: float) -> float:
    """
    Convert water mole fraction in wet gas to lb/mmscf of dry gas.

    Args:
        y_H2O: Water mole fraction in wet gas (mol/mol)

    Returns:
        Water content in lb per mmscf of dry gas
    """
    if y_H2O <= 0:
        return 0.0
    if y_H2O >= 1:
        raise ValueError("y_H2O must be less than 1")

    # At SC: 1 mmscf = 2635.88 lbmol, MW_water = 18.015
    CONV_FACTOR = 47476.0
    return y_H2O / (1.0 - y_H2O) * CONV_FACTOR


def ppm_to_molality(ppm_NaCl: float) -> float:
    """
    Convert ppm (mg/kg solution) NaCl to molality.

    Args:
        ppm_NaCl: Salt concentration in ppm (mg NaCl per kg solution)

    Returns:
        Molality in mol NaCl / kg water
    """
    if ppm_NaCl <= 0:
        return 0.0
    # ppm = mg/kg solution = g/1000kg solution
    # For dilute solutions: kg solution ≈ kg water
    # More precisely: wt% = ppm / 10000
    wt_pct = ppm_NaCl / 10000.0
    return wt_pct_to_molality(wt_pct)


def molality_to_ppm(molality: float) -> float:
    """
    Convert molality to ppm (mg/kg solution) NaCl.

    Args:
        molality: mol NaCl / kg water

    Returns:
        ppm NaCl (mg NaCl per kg solution)
    """
    if molality <= 0:
        return 0.0
    wt_pct = molality_to_wt_pct(molality)
    return wt_pct * 10000.0


def kelvin_to_fahrenheit(T_K: float) -> float:
    """Convert Kelvin to Fahrenheit."""
    return (T_K - 273.15) * 9/5 + 32


def pascal_to_psia(P_Pa: float) -> float:
    """Convert Pascal to psia."""
    return P_Pa / 6894.757


# =============================================================================
# Flexible Unit Input Helpers
# =============================================================================
def parse_temperature(value: float, unit: str) -> float:
    """
    Convert temperature to Kelvin.

    Args:
        value: Temperature value
        unit: Unit string ('K', 'C', 'F', 'degC', 'degF', 'celsius', 'fahrenheit', 'kelvin')

    Returns:
        Temperature in Kelvin
    """
    unit_lower = unit.lower().strip()
    if unit_lower in ('k', 'kelvin'):
        return value
    elif unit_lower in ('c', 'degc', 'celsius', '°c'):
        return celsius_to_kelvin(value)
    elif unit_lower in ('f', 'degf', 'fahrenheit', '°f'):
        return fahrenheit_to_kelvin(value)
    else:
        raise ValueError(f"Unknown temperature unit: {unit}. Use 'K', 'C', or 'F'")


def parse_pressure(value: float, unit: str) -> float:
    """
    Convert pressure to Pascal.

    Args:
        value: Pressure value
        unit: Unit string ('Pa', 'bar', 'bara', 'psia', 'psi', 'MPa', 'kPa')

    Returns:
        Pressure in Pascal
    """
    unit_lower = unit.lower().strip()
    if unit_lower in ('pa', 'pascal'):
        return value
    elif unit_lower in ('bar', 'bara'):
        return bar_to_pascal(value)
    elif unit_lower in ('psia', 'psi'):
        return psia_to_pascal(value)
    elif unit_lower == 'mpa':
        return value * 1e6
    elif unit_lower == 'kpa':
        return value * 1e3
    else:
        raise ValueError(f"Unknown pressure unit: {unit}. Use 'bar', 'psia', 'Pa', 'MPa', or 'kPa'")


def parse_salinity(value: float, unit: str) -> float:
    """
    Convert salinity to molality.

    Args:
        value: Salinity value
        unit: Unit string ('molal', 'mol/kg', 'ppm', 'wt%', 'wtpct', 'mg/L', 'g/L')

    Returns:
        Salinity in molality (mol NaCl / kg water)
    """
    unit_lower = unit.lower().strip()
    if unit_lower in ('molal', 'mol/kg', 'm'):
        return value
    elif unit_lower in ('ppm', 'mg/kg', 'mg/l'):
        return ppm_to_molality(value)
    elif unit_lower in ('wt%', 'wtpct', 'wt pct', 'weight%', 'percent'):
        return wt_pct_to_molality(value)
    elif unit_lower in ('g/l', 'g/kg'):
        # g/L ≈ g/kg for dilute solutions ≈ 1000 ppm
        return ppm_to_molality(value * 1000)
    else:
        raise ValueError(f"Unknown salinity unit: {unit}. Use 'molal', 'ppm', or 'wt%'")


# =============================================================================
# Binary VLE Calculator (Decoupled S&W Approach)
# =============================================================================
class SWBinaryVLE:
    """
    Soreide-Whitson VLE calculator for binary gas-water/brine systems.

    Uses DECOUPLED approach per Curtis Whitson:
    - Gas solubility: kij_AQ for BOTH phases
    - Water content: kij_NA for BOTH phases

    Framework modes:
    - 'proposed' (default): MC-3 alpha + proposed freshwater kij + Sechenov for all gases
    - 'sw_original': S&W alpha (with embedded salinity) + S&W kij + Eq 8 for all gases

    Example:
        vle = SWBinaryVLE('H2', salinity_molal=0.0)
        x_gas = vle.calc_gas_solubility(T_K, P_Pa)
        y_H2O = vle.calc_water_content(T_K, P_Pa)
    """

    def __init__(self, gas: str, salinity_molal: float = 0.0,
                 framework: str = 'proposed'):
        """
        Initialize VLE calculator.

        Args:
            gas: Gas name (e.g., 'H2', 'CO2', 'CH4')
            salinity_molal: Salt concentration in mol/kg water
            framework: 'proposed' (MC-3 alpha, freshwater kij, Sechenov routing),
                       'sw_original' (S&W alpha, embedded salinity kij, Eq 8), or
                       'dropin' (S&W alpha, new freshwater kij, embedded delta)
        """
        if gas not in COMPONENTS:
            raise ValueError(f"Unknown gas: {gas}. Supported: {list(COMPONENTS.keys())}")
        if framework not in ('proposed', 'sw_original', 'dropin'):
            raise ValueError(f"Unknown framework: {framework}. "
                             "Use 'proposed', 'sw_original', or 'dropin'")

        self.gas = gas
        self.salinity = salinity_molal
        self.framework = framework
        self.Tc = np.array([COMPONENTS['H2O'].Tc, COMPONENTS[gas].Tc])
        self.Pc = np.array([COMPONENTS['H2O'].Pc, COMPONENTS[gas].Pc])
        self.omega = np.array([COMPONENTS['H2O'].omega, COMPONENTS[gas].omega])
        self.sqrt2 = np.sqrt(2.0)
        self._last_vlle_warning = False

    def _calc_alpha(self, T_K: float) -> np.ndarray:
        """Calculate alpha for both components (framework-dependent for water)."""
        alpha = np.zeros(2)
        Tr_w = T_K / self.Tc[0]
        if self.framework in ('sw_original', 'dropin'):
            alpha[0] = alpha_water_soreide(Tr_w, self.salinity)
        else:
            alpha[0] = alpha_water_mc3(Tr_w)
        Tr_g = T_K / self.Tc[1]
        alpha[1] = alpha_standard_pr(Tr_g, self.omega[1])
        return alpha

    def _calc_eos_params(self, T_K: float, P_Pa: float) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate EOS A and B parameters."""
        alpha = self._calc_alpha(T_K)
        ai = OMEGA_A * (R_GAS * self.Tc)**2 * alpha / self.Pc
        bi = OMEGA_B * R_GAS * self.Tc / self.Pc
        Ai = ai * P_Pa / (R_GAS * T_K)**2
        Bi = bi * P_Pa / (R_GAS * T_K)
        return Ai, Bi

    def _calc_Aij_matrix(self, Ai: np.ndarray, kij: float) -> np.ndarray:
        """Calculate binary interaction matrix."""
        Aij = np.zeros((2, 2))
        Aij[0, 0] = Ai[0]
        Aij[1, 1] = Ai[1]
        Aij[0, 1] = (1.0 - kij) * np.sqrt(Ai[0] * Ai[1])
        Aij[1, 0] = Aij[0, 1]
        return Aij

    def _calc_K_with_kij(self, T_K: float, P_Pa: float, x: np.ndarray,
                         y: np.ndarray, kij: float) -> np.ndarray:
        """Calculate K-values using SAME kij for both phases."""
        Ai, Bi = self._calc_eos_params(T_K, P_Pa)
        Aij = self._calc_Aij_matrix(Ai, kij)

        # Liquid phase
        A_liq = np.einsum('i,j,ij', x, x, Aij)
        B_liq = np.dot(x, Bi)
        Z_liq = min(solve_cubic_eos(A_liq, B_liq))
        sum_xA = np.array([np.dot(x, Aij[i, :]) for i in range(2)])
        phi_liq = np.array([calc_fugacity_coeff(Z_liq, A_liq, B_liq,
                           Bi[i]/B_liq, sum_xA[i]/A_liq) for i in range(2)])

        # Vapor phase - SAME kij
        A_vap = np.einsum('i,j,ij', y, y, Aij)
        B_vap = np.dot(y, Bi)
        Z_vap = max(solve_cubic_eos(A_vap, B_vap))
        sum_yA = np.array([np.dot(y, Aij[i, :]) for i in range(2)])
        phi_vap = np.array([calc_fugacity_coeff(Z_vap, A_vap, B_vap,
                           Bi[i]/B_vap, sum_yA[i]/A_vap) for i in range(2)])

        return phi_liq / (phi_vap + 1e-30)

    def calc_gas_solubility(self, T_K: float, P_Pa: float, max_iter: int = 200,
                            salinity_method: str = 'auto') -> float:
        """
        Calculate gas mole fraction in aqueous phase using kij_AQ.

        Uses brentq root-finding on the fugacity balance (pure-gas vapor assumption).
        The pure-gas assumption (y_gas=1) is consistent with how kij_AQ was regressed
        and keeps the aqueous BIP loop independent of kij_NA.
        Returns np.nan when no two-phase VLE solution exists (e.g., above the mixture
        phase boundary for near-critical gases like H2S at high T).

        Args:
            T_K: Temperature in Kelvin
            P_Pa: Pressure in Pascal
            max_iter: Unused, retained for API compatibility
            salinity_method: 'auto' (framework default), 'sechenov', or 'embedded'

        Returns:
            x_gas: Gas mole fraction in aqueous phase (np.nan if no VLE solution)
        """
        fw = self.framework

        if salinity_method == 'auto':
            # Default salinity handling depends on framework
            if fw == 'proposed':
                # Proposed: freshwater kij + Sechenov for ALL gases when salinity > 0
                kij_aq = get_kij_aq(self.gas, T_K, 0.0, framework=fw)
                x_gas = self._calc_x_with_kij(T_K, P_Pa, kij_aq)
                if self.salinity > 0:
                    P_bar = P_Pa / 1e5
                    ks = get_sechenov_ks(self.gas, T_K, self.salinity, P_bar, framework=fw)
                    x_gas *= 10**(-ks * self.salinity)
            elif fw == 'dropin':
                # Drop-in: freshwater kij (dropin) + embedded delta for brine
                kij_aq = get_kij_aq(self.gas, T_K, 0.0, framework=fw)
                if self.salinity > 0 and self.gas in EMBEDDED_SALINITY_PARAMS_DROPIN:
                    delta = calc_embedded_delta_kij(
                        self.gas, T_K, self.salinity,
                        params=EMBEDDED_SALINITY_PARAMS_DROPIN[self.gas])
                    kij_aq += delta
                x_gas = self._calc_x_with_kij(T_K, P_Pa, kij_aq)
            else:
                # S&W original: salinity embedded in kij for HC/CO2/N2
                kij_aq = get_kij_aq(self.gas, T_K, self.salinity, framework=fw)
                x_gas = self._calc_x_with_kij(T_K, P_Pa, kij_aq)
                # H2 and H2S lack embedded salinity → apply Sechenov
                if self.salinity > 0 and self.gas not in _SW_GASES_WITH_EMBEDDED_SALINITY:
                    P_bar = P_Pa / 1e5
                    ks = get_sechenov_ks(self.gas, T_K, self.salinity, P_bar, framework=fw)
                    x_gas *= 10**(-ks * self.salinity)

        elif salinity_method == 'embedded':
            # Embedded BIP: kij_fw(T) + delta(T,m)
            kij_fw = get_kij_aq(self.gas, T_K, 0.0, framework=fw)
            if fw == 'sw_original' and self.gas in _SW_GASES_WITH_EMBEDDED_SALINITY:
                # S&W original embedded gases: use S&W kij directly (salinity already in kij)
                kij_aq = get_kij_aq(self.gas, T_K, self.salinity, framework=fw)
            elif fw == 'dropin' and self.gas in EMBEDDED_SALINITY_PARAMS_DROPIN:
                delta = calc_embedded_delta_kij(
                    self.gas, T_K, self.salinity,
                    params=EMBEDDED_SALINITY_PARAMS_DROPIN[self.gas])
                kij_aq = kij_fw + delta
            else:
                # Use Paper 2 embedded delta (proposed framework)
                delta = calc_embedded_delta_kij(self.gas, T_K, self.salinity)
                kij_aq = kij_fw + delta
            x_gas = self._calc_x_with_kij(T_K, P_Pa, kij_aq)

        elif salinity_method == 'sechenov':
            # Explicit Sechenov: freshwater kij + post-solve correction
            kij_aq = get_kij_aq(self.gas, T_K, 0.0, framework=fw)
            x_gas = self._calc_x_with_kij(T_K, P_Pa, kij_aq)
            if self.salinity > 0:
                P_bar = P_Pa / 1e5
                ks = get_sechenov_ks(self.gas, T_K, self.salinity, P_bar, framework=fw)
                x_gas *= 10**(-ks * self.salinity)
        else:
            raise ValueError(f"Unknown salinity_method: {salinity_method}")

        return x_gas

    def calc_water_content(self, T_K: float, P_Pa: float, max_iter: int = 200) -> float:
        """
        Calculate water mole fraction in gas phase using kij_NA.

        Args:
            T_K: Temperature in Kelvin
            P_Pa: Pressure in Pascal
            max_iter: Maximum iterations

        Returns:
            y_H2O: Water mole fraction in gas phase
        """
        kij_na = get_kij_na(self.gas, T_K)
        return self.calc_water_content_with_kij(T_K, P_Pa, kij_na, max_iter)

    def calc_water_content_with_kij(self, T_K: float, P_Pa: float, kij_na: float,
                                     max_iter: int = 200) -> float:
        """
        Calculate water mole fraction in gas phase using specified kij_NA.

        Args:
            T_K: Temperature in Kelvin
            P_Pa: Pressure in Pascals
            kij_na: Binary interaction parameter for non-aqueous phase
            max_iter: Maximum iterations

        Returns:
            y_H2O: Water mole fraction in gas phase
        """
        x = np.array([0.999, 0.001])
        y = np.array([0.02, 0.98])

        for iteration in range(max_iter):
            x = np.clip(x, 1e-14, 1.0 - 1e-14)
            x = x / np.sum(x)
            y = np.clip(y, 1e-14, 1.0 - 1e-14)
            y = y / np.sum(y)

            K = self._calc_K_with_kij(T_K, P_Pa, x, y, kij_na)

            y_new = K * x
            y_new = np.clip(y_new, 1e-14, 1.0 - 1e-14)
            y_new = y_new / np.sum(y_new)

            error = np.max(np.abs(y_new - y))
            if error < 1e-10:
                return y[0]

            damp = 0.4
            y = y + damp * (y_new - y)
            y = y / np.sum(y)

            x_new = y / (K + 1e-30)
            x_new = np.clip(x_new, 1e-14, 1.0 - 1e-14)
            x_new = x_new / np.sum(x_new)
            if x_new[0] < 0.5:
                x_new = np.array([0.98, 0.02])
            x = x + damp * (x_new - x)
            x = x / np.sum(x)

        return y[0]

    # Aliases for backward compatibility
    def calc_x_gas(self, T_K: float, P_Pa: float, kij: Optional[float] = None) -> float:
        """Alias for calc_gas_solubility. If kij provided, uses that value."""
        if kij is not None:
            return self._calc_x_with_kij(T_K, P_Pa, kij)
        return self.calc_gas_solubility(T_K, P_Pa)

    def calc_x_H2(self, T_K: float, P_Pa: float, kij: float) -> float:
        """Calculate H2 solubility with specified kij (for regression)."""
        return self._calc_x_with_kij(T_K, P_Pa, kij)

    def _calc_x_with_kij(self, T_K: float, P_Pa: float, kij: float) -> float:
        """
        Calculate gas solubility with specified kij_AQ using brentq root-finding.

        Solves for x_gas where fugacity of gas in liquid equals fugacity in vapor.
        Uses pure gas vapor phase (y_gas = 1.0) assumption — consistent with how
        kij_AQ was regressed and keeps the aqueous BIP independent of kij_NA.
        Scans for first valid bracket to handle high-solubility and near-critical cases.
        Returns np.nan if no root can be found.
        """
        from scipy.optimize import brentq

        # Pure gas vapor phase
        y = np.array([0.0, 1.0])
        phi_V = self._calc_fugacity_coeff_with_kij(T_K, P_Pa, y, kij, phase='vapor')
        f_gas_V = y[1] * phi_V[1] * P_Pa

        def objective(x_gas):
            x = np.array([1.0 - x_gas, x_gas])
            phi_L = self._calc_fugacity_coeff_with_kij(T_K, P_Pa, x, kij, phase='liquid')
            f_gas_L = x[1] * phi_L[1] * P_Pa
            return f_gas_L - f_gas_V

        # Try standard narrow bracket first (covers most cases efficiently)
        try:
            return brentq(objective, 1e-10, 0.2)
        except ValueError:
            pass

        # Bracket scan: probe objective at logarithmically-spaced x values
        # to find the first sign change (physical root at lowest x_gas)
        scan_pts = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 5e-3,
                    0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40,
                    0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
        vals = []
        for xp in scan_pts:
            try:
                vals.append(objective(xp))
            except Exception:
                vals.append(np.nan)

        for i in range(len(scan_pts) - 1):
            if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]):
                if vals[i] * vals[i + 1] < 0:
                    try:
                        return brentq(objective, scan_pts[i], scan_pts[i + 1])
                    except Exception:
                        continue

        return np.nan

    def _calc_fugacity_coeff_with_kij(self, T_K: float, P_Pa: float,
                                       x: np.ndarray, kij: float,
                                       phase: str = 'liquid') -> np.ndarray:
        """Calculate fugacity coefficients for both components with specified kij."""
        Ai, Bi = self._calc_eos_params(T_K, P_Pa)
        Aij = self._calc_Aij_matrix(Ai, kij)

        A_mix = sum(x[i] * x[j] * Aij[i, j] for i in range(2) for j in range(2))
        B_mix = sum(x[i] * Bi[i] for i in range(2))

        roots = solve_cubic_eos(A_mix, B_mix)
        Z = roots[0] if phase == 'liquid' else roots[-1]

        phi = np.zeros(2)
        for i in range(2):
            sum_xA_i = sum(x[j] * Aij[i, j] for j in range(2))
            phi[i] = calc_fugacity_coeff(Z, A_mix, B_mix, Bi[i]/B_mix if B_mix > 0 else 0,
                                          sum_xA_i/A_mix if A_mix > 0 else 0)
        return phi

    def calc_y_H2O(self, T_K: float, P_Pa: float, kij_na: float) -> float:
        """Alias for calc_water_content_with_kij."""
        return self.calc_water_content_with_kij(T_K, P_Pa, kij_na)


# =============================================================================
# Gas-Gas BIP Database (Literature Values)
# =============================================================================
# Sources: GPSA Engineering Data Book, Knapp et al., various EOS studies
# Convention: GAS_GAS_BIPS[(gas_a, gas_b)] = kij
# Unspecified pairs default to 0.0

GAS_GAS_BIPS = {
    ('CH4', 'CO2'): 0.12,      ('CO2', 'C2H6'): 0.13,
    ('CO2', 'C3H8'): 0.135,    ('CO2', 'N2'): -0.02,
    ('CO2', 'H2S'): 0.097,     ('CO2', 'H2'): 0.0,
    ('CO2', 'nC4H10'): 0.13,   ('CO2', 'iC4H10'): 0.13,
    ('CH4', 'C2H6'): 0.0026,   ('CH4', 'C3H8'): 0.014,
    ('CH4', 'N2'): 0.036,      ('CH4', 'H2S'): 0.08,
    ('CH4', 'H2'): 0.0,        ('CH4', 'nC4H10'): 0.02,
    ('CH4', 'iC4H10'): 0.02,   ('H2', 'N2'): 0.0,
    ('H2S', 'N2'): 0.17,       ('C2H6', 'N2'): 0.04,
    ('C3H8', 'N2'): 0.08,      ('H2', 'H2S'): 0.0,
    ('C2H6', 'H2S'): 0.085,    ('C3H8', 'H2S'): 0.08,
    ('C2H6', 'H2'): 0.0,       ('C3H8', 'H2'): 0.0,
    ('H2', 'nC4H10'): 0.0,     ('H2', 'iC4H10'): 0.0,
    ('C2H6', 'C3H8'): 0.001,   ('C2H6', 'nC4H10'): 0.01,
    ('C3H8', 'nC4H10'): 0.003,
}


def get_gas_gas_bip(gas_a: str, gas_b: str) -> float:
    """Get gas-gas BIP from database. Returns 0.0 for unknown pairs."""
    if gas_a == gas_b:
        return 0.0
    for key in [(gas_a, gas_b), (gas_b, gas_a)]:
        if key in GAS_GAS_BIPS:
            return GAS_GAS_BIPS[key]
    return 0.0


# =============================================================================
# Multi-Component Flash Calculator (Curtis Whitson scheme, Feb 2026)
# =============================================================================
# =============================================================================
# S&W-Specific K-Value Initialization
# =============================================================================
# Component-specific K-value parameters fitted to S&W binary VLE results.
# All light gases use the Cross form:
#   ln(K) = a + b(Tc/T) + c·ln(P/Pc) + d(Tc/T)² + e·ln(P/Pc)² + f(Tc/T)·ln(P/Pc)
# Heavy HCs (C5+) use LogLinear + floor(K=10):
#   ln(K) = a + b(Tc/T) + c·ln(P/Pc), K = max(K, 10)
# Water uses a universal 6-parameter T-P correlation.
#
# Fitted across: 14.7–15000 psia (1–1034 bar), 32–300°F (273–422 K)
# using 500 binary VLE evaluations per component (20 T × 25 P, log-spaced P).
# Fitting: robust Huber loss (f_scale=0.5) with 10 random restarts.

_SW_KVALUE_PARAMS = {
    # Light gases: Cross form [a, b, c, d, e, f], MARE on K
    # Fitted using proposed Paper 2 freshwater kij_AQ forms (Feb 2026)
    'H2':     ([6.4295, 25.5844, -0.5985, -50.0000, -0.0007, -3.2598], 11.7),
    'CO2':    ([-5.9974, 26.3804, -0.9380, -16.3941, 0.0607, 0.3688], 13.7),
    'N2':     ([1.4998, 36.4929, -0.6419, -50.0000, 0.0110, -0.6328], 8.4),
    'H2S':    ([-1.0549, 7.7334, -1.3646, -3.4035, 0.0850, 0.8651], 18.8),
    'CH4':    ([-2.7107, 36.5276, -0.7040, -33.2023, 0.0312, -0.1419], 8.5),
    'C2H6':   ([-9.8804, 40.0208, -0.9108, -22.8267, 0.0756, 0.4366], 10.4),
    'C3H8':   ([-9.8805, 33.3750, -1.0803, -15.1332, 0.0836, 0.6959], 12.8),
    'iC4H10': ([-8.8362, 29.1657, -1.0755, -11.0413, 0.0685, 0.7280], 17.9),
    'nC4H10': ([-8.4159, 27.0161, -1.0886, -10.2525, 0.0613, 0.7515], 21.6),
}

_SW_KVALUE_HEAVY = {
    # Heavy HCs: LogLinear form [a, b, c] + floor=10
    'iC5H12':  [6.5325, 3.4733, -0.0386],
    'nC5H12':  [6.2726, 3.7905, -0.0083],
    'nC6H14':  [4.5498, 6.4018, 0.1051],
    'nC7H16':  [2.6475, 8.9734, 0.1954],
    'nC8H18':  [0.3846, 11.9724, 0.2711],
    'nC10H22': [3.8490, -1.1163, -0.0764],
}

# Water K-value: ln(K_H2O) = a + b/T + c·ln(P) + d/T² + e·ln(P)/T + f·ln(P)²
_SW_KVALUE_WATER = [34.3692, -1873.3949, -3.2671, -555202.1539, 42.1759, 0.0804]


def _sw_kvalue_init(names: list, Tc: np.ndarray, Pc: np.ndarray,
                    omega: np.ndarray, T_K: float, P_Pa: float) -> np.ndarray:
    """
    S&W-specific K-value initialization for gas-water flash.

    Returns array of K-values for all components, using correlations fitted
    to binary S&W VLE results. Eliminates the K < 1 violations that plague
    Wilson initialization for gas-water systems.

    Args:
        names: component name list (must include 'H2O')
        Tc, Pc, omega: arrays of critical properties
        T_K: temperature (K)
        P_Pa: pressure (Pa)

    Returns:
        K: array of K-values (one per component)
    """
    nc = len(names)
    K = np.ones(nc)
    T = T_K
    P = P_Pa

    for i, name in enumerate(names):
        if name == 'H2O':
            # Universal water K-value (6-param T-P fit)
            a, b, c, d, e, f = _SW_KVALUE_WATER
            K[i] = np.exp(a + b/T + c*np.log(P) + d/T**2
                          + e*np.log(P)/T + f*np.log(P)**2)

        elif name in _SW_KVALUE_PARAMS:
            # Light gas: Cross form
            p, _ = _SW_KVALUE_PARAMS[name]
            Tr_inv = Tc[i] / T
            lnPr = np.log(P / Pc[i])
            K[i] = np.exp(p[0] + p[1]*Tr_inv + p[2]*lnPr
                          + p[3]*Tr_inv**2 + p[4]*lnPr**2
                          + p[5]*Tr_inv*lnPr)

        elif name in _SW_KVALUE_HEAVY:
            # Heavy HC: LogLinear + floor
            p = _SW_KVALUE_HEAVY[name]
            Tr_inv = Tc[i] / T
            lnPr = np.log(P / Pc[i])
            K[i] = max(np.exp(p[0] + p[1]*Tr_inv + p[2]*lnPr), 10.0)

        else:
            # Fallback: standard Wilson for unknown components
            K[i] = (Pc[i] / P) * np.exp(
                5.373 * (1.0 + omega[i]) * (1.0 - Tc[i] / T))

    return np.clip(K, 1e-10, 1e10)


class SWMultiComponentFlash:
    """
    Full N-component flash per Curtis Whitson's confirmed S&W scheme:

    Flash 1: All gas-water BIPs = kij_AQ, standard Rachford-Rice → take AQUEOUS phase
    Flash 2: All gas-water BIPs = kij_NA, standard Rachford-Rice → take NON-AQUEOUS phase
    True K-values: K_i = y_i(Flash 2) / x_i(Flash 1)

    Gas-gas BIPs from literature used in both flashes.

    Usage:
        flash = SWMultiComponentFlash(['H2O', 'CH4', 'CO2', 'N2'])
        result = flash.calc_equilibrium(T_K, P_Pa, z)
    """

    def __init__(self, component_names: List[str], salinity_molal: float = 0.0,
                 framework: str = 'proposed'):
        if 'H2O' not in component_names:
            raise ValueError("H2O must be in component list")
        if framework not in ('proposed', 'sw_original', 'dropin'):
            raise ValueError(f"Unknown framework: {framework}. "
                             "Use 'proposed', 'sw_original', or 'dropin'")
        self.names = list(component_names)
        self.nc = len(self.names)
        self.salinity = salinity_molal
        self.framework = framework
        self.iw = self.names.index('H2O')
        self.Tc = np.array([COMPONENTS[n].Tc for n in self.names])
        self.Pc = np.array([COMPONENTS[n].Pc for n in self.names])
        self.omega = np.array([COMPONENTS[n].omega for n in self.names])
        self.Tb = np.array([COMPONENTS[n].Tb for n in self.names])
        self._last_vlle_warning = False

    def _calc_alpha(self, T_K: float) -> np.ndarray:
        alpha = np.zeros(self.nc)
        for i in range(self.nc):
            if self.names[i] == 'H2O':
                Tr_w = T_K / self.Tc[i]
                if self.framework in ('sw_original', 'dropin'):
                    alpha[i] = alpha_water_soreide(Tr_w, self.salinity)
                else:
                    alpha[i] = alpha_water_mc3(Tr_w)
            else:
                alpha[i] = alpha_standard_pr(T_K / self.Tc[i], self.omega[i])
        return alpha

    def _calc_ai_bi(self, T_K: float) -> Tuple[np.ndarray, np.ndarray]:
        alpha = self._calc_alpha(T_K)
        ai = OMEGA_A * (R_GAS * self.Tc)**2 * alpha / self.Pc
        bi = OMEGA_B * R_GAS * self.Tc / self.Pc
        return ai, bi

    def build_kij_matrix(self, T_K: float, mode: str = 'AQ') -> np.ndarray:
        """Build N×N kij matrix. mode='AQ' or 'NA' for gas-water pairs."""
        kij = np.zeros((self.nc, self.nc))
        for i in range(self.nc):
            for j in range(i + 1, self.nc):
                ni, nj = self.names[i], self.names[j]
                if ni == 'H2O' or nj == 'H2O':
                    gas = nj if ni == 'H2O' else ni
                    if mode == 'AQ':
                        val = get_kij_aq(gas, T_K, self.salinity,
                                         framework=self.framework)
                    else:
                        val = get_kij_na(gas, T_K)
                else:
                    val = get_gas_gas_bip(ni, nj)
                kij[i, j] = val
                kij[j, i] = val
        return kij

    def calc_fugacity_coefficients(self, T_K: float, P_Pa: float, comp: np.ndarray,
                                    kij_matrix: np.ndarray, phase: str = 'liquid') -> np.ndarray:
        """Calculate fugacity coefficients for all components."""
        ai, bi = self._calc_ai_bi(T_K)
        RT = R_GAS * T_K
        Ai = ai * P_Pa / RT**2
        Bi = bi * P_Pa / RT

        # Vectorized mixing rule (replaces O(nc^2) Python loop)
        sqrt_Ai = np.sqrt(Ai)
        Aij = np.outer(sqrt_Ai, sqrt_Ai) * (1.0 - kij_matrix)

        A_mix = np.einsum('i,j,ij', comp, comp, Aij)
        B_mix = np.dot(comp, Bi)

        if B_mix < 1e-15 or A_mix < 1e-15:
            return np.ones(self.nc)

        roots = solve_cubic_eos(A_mix, B_mix)
        Z = roots[0] if phase == 'liquid' else roots[-1]

        # Vectorized fugacity coefficient calculation
        sum_xA = Aij @ comp  # matrix-vector product replaces per-component loop
        Bi_over_B = Bi / B_mix
        sum_xA_over_A = sum_xA / A_mix

        phi = np.zeros(self.nc)
        for i in range(self.nc):
            phi[i] = calc_fugacity_coeff(Z, A_mix, B_mix, Bi_over_B[i], sum_xA_over_A[i])
        return phi

    def _calc_fugacity_fast(self, comp: np.ndarray, Ai: np.ndarray, Bi: np.ndarray,
                            sqrt_Ai: np.ndarray, onemk: np.ndarray,
                            phase: str = 'liquid') -> np.ndarray:
        """
        Fast fugacity coefficient calculation with precomputed T/P-dependent quantities.

        Used by flash_tp and _calc_x_with_kij to avoid recomputing ai, bi, Ai, Bi
        on every iteration when only compositions change.
        """
        Aij = np.outer(sqrt_Ai, sqrt_Ai) * onemk

        A_mix = np.einsum('i,j,ij', comp, comp, Aij)
        B_mix = np.dot(comp, Bi)

        if B_mix < 1e-15 or A_mix < 1e-15:
            return np.ones(len(comp))

        roots = solve_cubic_eos(A_mix, B_mix)
        Z = roots[0] if phase == 'liquid' else roots[-1]

        sum_xA = Aij @ comp
        Bi_over_B = Bi / B_mix
        sum_xA_over_A = sum_xA / A_mix

        phi = np.zeros(len(comp))
        for i in range(len(comp)):
            phi[i] = calc_fugacity_coeff(Z, A_mix, B_mix, Bi_over_B[i], sum_xA_over_A[i])
        return phi

    def _wilson_k_values(self, T_K: float, P_Pa: float) -> np.ndarray:
        """
        S&W-specific K-value initialization for gas-water flash.

        Replaces the standard Wilson correlation with component-specific fits
        derived from binary S&W VLE results across 1–1034 bar, 273–422 K.

        The standard Wilson correlation K = (Pc/P)·exp(5.373(1+ω)(1-Tc/T))
        was developed for hydrocarbon mixtures and gives qualitatively wrong
        K-values (K < 1) for many gas-water systems at moderate conditions.
        For example, Wilson predicts K_H2S < 1 at 44% of tested conditions,
        even though H2S always prefers the vapor phase in gas-water equilibria.

        These replacement correlations ensure K_gas > 1 for all gas components
        at all conditions in the valid range, providing a reliable starting
        point for successive substitution convergence.

        Light gases (H2 through nC4H10) use the 'Cross' form:
          ln(K) = a + b(Tc/T) + c·ln(P/Pc) + d(Tc/T)² + e·ln(P/Pc)² + f(Tc/T)·ln(P/Pc)
        where the cross term f(Tc/T)·ln(P/Pc) captures the T-P interaction
        characteristic of gas-water equilibria.

        Heavy hydrocarbons (C5+) use a simpler LogLinear form with a floor:
          ln(K) = a + b(Tc/T) + c·ln(P/Pc),  K = max(K, 10)
        These always have K >> 1, so quantitative accuracy matters less.

        Water uses a universal 6-parameter fit in T and P.

        See: code/fit_flash_kvalues.py for the reproducible fitting procedure.
        See: data/flash_kvalue_init_method.md for methodology documentation.
        """
        return _sw_kvalue_init(self.names, self.Tc, self.Pc, self.omega,
                               T_K, P_Pa)

    def calc_gamma(self, T_K: float, P_bar: float = 100.0,
                   salinity_molal: Optional[float] = None,
                   ks_override: Optional[Dict[str, float]] = None) -> np.ndarray:
        """
        Calculate Sechenov activity coefficients for gamma-phi flash.

        γ_i = 10^(ks_i × m) for gases
        γ_H2O = 1.0

        Uses get_sechenov_ks() for gas-specific ks routing:
        - CO2: Duan & Sun 2003 Pitzer model
        - H2S: Akinfiev et al. 2016 Pitzer model
        - All others (HCs, N2, H2): S&W Equation 8 with Tb

        Args:
            T_K: Temperature in Kelvin
            P_bar: Pressure in bar (needed for Duan 2003 CO2). Default 100.
            salinity_molal: Override salinity (defaults to self.salinity)
            ks_override: Optional dict of {gas_name: ks_value} to override
                         the default routing for specific gases

        Returns:
            Array of activity coefficients (one per component)
        """
        m = salinity_molal if salinity_molal is not None else self.salinity
        gamma = np.ones(self.nc)

        if m <= 0:
            return gamma

        for i in range(self.nc):
            if self.names[i] != 'H2O':
                if ks_override and self.names[i] in ks_override:
                    ks = ks_override[self.names[i]]
                else:
                    ks = get_sechenov_ks(self.names[i], T_K, m, P_bar,
                                        framework=self.framework)
                gamma[i] = 10.0**(ks * m)

        return gamma

    def flash_tp(self, T_K: float, P_Pa: float, z: np.ndarray,
                 mode: str = 'AQ', gamma: Optional[np.ndarray] = None,
                 max_iter: int = 200,
                 tol: float = 1e-10) -> Tuple[float, np.ndarray, np.ndarray, bool]:
        """
        Full TP flash using successive substitution with robust RR solver.

        Args:
            T_K: Temperature in Kelvin
            P_Pa: Pressure in Pascal
            z: Feed composition
            mode: 'AQ' for aqueous phase BIPs, 'NA' for non-aqueous
            gamma: Optional activity coefficient array for gamma-phi flash.
                   If provided, K_i = γ_i × φ_i^L / φ_i^V.
                   Only meaningful for mode='AQ' (salt doesn't enter vapor).
            max_iter: Maximum SS iterations
            tol: Convergence tolerance on K-values

        Returns: V (vapor fraction), x (liquid), y (vapor), converged (bool)
        """
        z = np.asarray(z, dtype=float)
        z = z / np.sum(z)
        kij_matrix = self.build_kij_matrix(T_K, mode)

        # Precompute T/P-dependent quantities (constant across SS iterations)
        ai, bi = self._calc_ai_bi(T_K)
        RT = R_GAS * T_K
        Ai = ai * P_Pa / RT**2
        Bi = bi * P_Pa / RT
        sqrt_Ai = np.sqrt(Ai)
        onemk = 1.0 - kij_matrix

        # Activity coefficients (default = 1 = no salinity effect)
        if gamma is None:
            gamma_eff = np.ones(self.nc)
        else:
            gamma_eff = np.asarray(gamma, dtype=float)

        # Initialize K-values (Wilson + gamma for initial estimate)
        K = self._wilson_k_values(T_K, P_Pa) * gamma_eff
        K[self.iw] = min(K[self.iw], 0.01)

        converged = False
        for it in range(max_iter):
            # Robust RR solver (Nielsen & Lia 2022)
            V, x, y = solve_rachford_rice(z, K)
            x = np.clip(x, 1e-15, None)
            y = np.clip(y, 1e-15, None)
            x = x / np.sum(x)
            y = y / np.sum(y)

            phi_L = self._calc_fugacity_fast(x, Ai, Bi, sqrt_Ai, onemk, 'liquid')
            phi_V = self._calc_fugacity_fast(y, Ai, Bi, sqrt_Ai, onemk, 'vapor')

            # Gamma-phi K-value: K_i = γ_i × φ_i^L / φ_i^V
            K_new = np.clip(gamma_eff * phi_L / (phi_V + 1e-30), 1e-10, 1e10)

            if np.max(np.abs(K_new / K - 1.0)) < tol:
                converged = True
                K = K_new
                break

            damp = 0.7 if it < 20 else 0.9
            K = K * (K_new / K)**damp

        # Final compositions with converged K
        V, x, y = solve_rachford_rice(z, K)
        x = np.clip(x, 1e-15, None)
        y = np.clip(y, 1e-15, None)
        x = x / np.sum(x)
        y = y / np.sum(y)
        return V, x, y, converged

    def check_vlle(self, T_K: float, P_Pa: float, x_ref: np.ndarray,
                    kij_matrix: np.ndarray) -> bool:
        """
        TPD-based VLLE screening (opt-in diagnostic).

        Tests pure-component trial compositions against a reference liquid phase.
        If any TPD < -0.01, flags potential additional liquid phase.

        Non-invasive: does not modify any computed values.

        Args:
            T_K: Temperature in Kelvin
            P_Pa: Pressure in Pascal
            x_ref: Reference liquid phase composition (from converged flash)
            kij_matrix: BIP matrix used for the flash

        Returns:
            True if VLLE is suspected, False otherwise
        """
        ai, bi = self._calc_ai_bi(T_K)
        RT = R_GAS * T_K
        Ai = ai * P_Pa / RT**2
        Bi = bi * P_Pa / RT
        sqrt_Ai = np.sqrt(Ai)
        onemk = 1.0 - kij_matrix

        phi_ref = self._calc_fugacity_fast(x_ref, Ai, Bi, sqrt_Ai, onemk, 'liquid')
        ln_x_ref = np.log(np.clip(x_ref, 1e-30, None))
        ln_phi_ref = np.log(np.clip(phi_ref, 1e-30, None))

        nc = len(x_ref)
        for k in range(nc):
            w = np.full(nc, 1e-10)
            w[k] = 1.0
            w = w / np.sum(w)

            phi_trial = self._calc_fugacity_fast(w, Ai, Bi, sqrt_Ai, onemk, 'liquid')
            ln_phi_trial = np.log(np.clip(phi_trial, 1e-30, None))

            tpd_k = (np.log(w[k]) + ln_phi_trial[k]) - (ln_x_ref[k] + ln_phi_ref[k])
            if tpd_k < -0.01:
                return True

        return False

    def calc_equilibrium(self, T_K: float, P_Pa: float, z: np.ndarray,
                         salinity_method: str = 'gamma_phi',
                         ks_override: Optional[Dict[str, float]] = None,
                         salinity_for_sechenov: float = 0.0) -> Dict:
        """
        Full S&W equilibrium per Curtis's scheme.

        Flash 1 (kij_AQ) → aqueous phase x
        Flash 2 (kij_NA) → non-aqueous phase y

        Salinity methods:
        - 'gamma_phi' (RECOMMENDED): Activity coefficient γ_i = 10^(ks×m)
          included inside the aqueous flash K-value update. Standard
          electrolyte thermodynamics. Exact mass balance.
        - 'explicit': Freshwater flash first, then post-solve Sechenov
          correction. Slight mass balance error for concentrated systems.
        - 'embedded': Use kij(csw) with embedded salinity (original S&W).
          Only available for HC/CO2/N2. H2 falls back to gamma_phi.

        Args:
            T_K: Temperature in Kelvin
            P_Pa: Pressure in Pascal
            z: Feed composition
            salinity_method: 'gamma_phi', 'explicit', or 'embedded'
            ks_override: Optional dict {gas_name: ks_value} for gases
                         where Eq 8 may not apply
            salinity_for_sechenov: (Deprecated) For backward compatibility,
                         if > 0 and salinity_method not specified, uses
                         explicit Sechenov.

        Returns dict with:
            'x_aq': Aqueous phase mole fractions (from AQ flash)
            'y_na': Non-aqueous phase mole fractions (from NA flash)
            'K_true': True K-values y_na/x_aq
            'V_aq', 'V_na': Vapor fractions from each flash
            'converged_aq', 'converged_na': Convergence flags
            'gamma': Activity coefficients used (if gamma_phi)
            'component_names': Component name list
        """
        z = np.asarray(z, dtype=float) / np.sum(z)

        # Handle backward compatibility
        if salinity_for_sechenov > 0 and salinity_method == 'gamma_phi':
            salinity_method = 'explicit'

        # Calculate gamma for gamma-phi method
        gamma_aq = None
        if salinity_method == 'gamma_phi' and self.salinity > 0:
            P_bar = P_Pa / 1e5
            gamma_aq = self.calc_gamma(T_K, P_bar=P_bar, ks_override=ks_override)

        # Flash 1: Aqueous phase (kij_AQ)
        V_aq, x_aq, y_aq, conv_aq = self.flash_tp(
            T_K, P_Pa, z, mode='AQ', gamma=gamma_aq)

        # Flash 2: Non-aqueous phase (kij_NA) — no gamma (salt stays in liquid)
        V_na, x_na, y_na, conv_na = self.flash_tp(T_K, P_Pa, z, mode='NA')

        K_true = np.where(x_aq > 1e-15, y_na / x_aq, 1e10)

        result = {
            'x_aq': x_aq, 'y_na': y_na, 'K_true': K_true,
            'V_aq': V_aq, 'V_na': V_na,
            'converged_aq': conv_aq, 'converged_na': conv_na,
            'component_names': self.names,
            'salinity_method': salinity_method,
            'vlle_warning': self._last_vlle_warning,
        }

        if gamma_aq is not None:
            result['gamma'] = gamma_aq

        # Explicit Sechenov post-correction (backward compat or explicit request)
        if salinity_method == 'explicit':
            m = salinity_for_sechenov if salinity_for_sechenov > 0 else self.salinity
            if m > 0:
                T_C = T_K - 273.15
                x_brine = x_aq.copy()
                for i in range(self.nc):
                    if self.names[i] != 'H2O':
                        if ks_override and self.names[i] in ks_override:
                            ks = ks_override[self.names[i]]
                        else:
                            ks = sw_equation_8_ks(T_C, self.Tb[i])
                        x_brine[i] = x_aq[i] * 10**(-ks * m)
                x_brine[self.iw] = 1.0 - np.sum(x_brine) + x_brine[self.iw]
                result['x_aq_brine'] = x_brine

        return result


# =============================================================================
# Main Multi-Component API
# =============================================================================
def calc_gas_brine_equilibrium(
    salinity_wt_pct: float,
    temperature_F: float,
    pressure_psia: float,
    y_CH4: float = 0.0,
    y_C2H6: float = 0.0,
    y_C3H8: float = 0.0,
    y_nC4H10: float = 0.0,
    y_iC4H10: float = 0.0,
    y_CO2: float = 0.0,
    y_N2: float = 0.0,
    y_H2S: float = 0.0,
    y_H2: float = 0.0,
    method: str = 'henry',
    salinity_method: str = 'gamma_phi',
    framework: str = 'proposed',
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    Calculate gas-brine equilibrium using Soreide-Whitson framework.

    METHODS (flash structure):
    - 'henry' (default): Binary decomposition + Henry's law.
      x_i = x_i_pure * y_i. Valid for dilute solutions (<5 mol%).

    - 'flash': Full multi-component Rachford-Rice flash per Curtis Whitson
      scheme. Flash 1 (kij_AQ) → aqueous phase. Flash 2 (kij_NA) →
      non-aqueous phase. More accurate for concentrated systems.

    SALINITY_METHOD (salting-out approach, used with method='flash'):
    - 'gamma_phi' (default, RECOMMENDED): Freshwater EOS + Sechenov activity
      coefficient γ_i = 10^(ks_i × m) inside K-value update. Standard
      electrolyte thermodynamics. Exact mass balance. Trivial addition to
      the already customised S&W flash.

    - 'explicit': Freshwater flash + post-solve Sechenov correction (Eq 8).
      Algebraically equivalent to gamma_phi for dilute systems.

    - 'embedded': Original S&W salinity in kij (Eqs 12-15). Available for
      HC/CO2/N2. H2 falls back to gamma_phi.

    Args:
        salinity_wt_pct: Brine salinity in weight percent NaCl
        temperature_F: Temperature in degrees Fahrenheit
        pressure_psia: Pressure in psia
        y_*: Mole fractions of each gas in dry gas (will be normalized)
        method: 'henry' (binary decomposition) or 'flash' (multi-component)
        salinity_method: 'gamma_phi', 'explicit', or 'embedded'
        framework: 'proposed', 'sw_original', or 'dropin'

    Returns:
        Tuple of:
        - x_gas: Dict of dissolved gas mole fractions in brine
        - water_content: Dict with 'y_H2O', 'stb_mmscf', 'lb_mmscf'
    """
    # Input validation
    if salinity_wt_pct < 0 or salinity_wt_pct >= 100:
        raise ValueError("Salinity must be between 0 and 100 wt%")
    if pressure_psia <= 0:
        raise ValueError("Pressure must be positive")

    # Convert units
    T_K = fahrenheit_to_kelvin(temperature_F)
    P_Pa = psia_to_pascal(pressure_psia)
    salinity_molal = wt_pct_to_molality(salinity_wt_pct)

    # Build gas composition dict
    gas_comp = {}
    if y_CH4 > 0: gas_comp['CH4'] = y_CH4
    if y_C2H6 > 0: gas_comp['C2H6'] = y_C2H6
    if y_C3H8 > 0: gas_comp['C3H8'] = y_C3H8
    if y_nC4H10 > 0: gas_comp['nC4H10'] = y_nC4H10
    if y_iC4H10 > 0: gas_comp['iC4H10'] = y_iC4H10
    if y_CO2 > 0: gas_comp['CO2'] = y_CO2
    if y_N2 > 0: gas_comp['N2'] = y_N2
    if y_H2S > 0: gas_comp['H2S'] = y_H2S
    if y_H2 > 0: gas_comp['H2'] = y_H2

    if not gas_comp:
        raise ValueError("At least one gas component must have non-zero mole fraction")

    # Normalize
    total = sum(gas_comp.values())
    gas_comp = {k: v/total for k, v in gas_comp.items()}

    # Calculate solubility for each gas
    if method == 'flash':
        # Full multi-component flash (Curtis Whitson scheme)
        comp_names = ['H2O'] + list(gas_comp.keys())
        z = np.zeros(len(comp_names))
        z[0] = 0.95  # 95% water feed
        gas_frac = 0.05
        for i, (gas, y_i) in enumerate(gas_comp.items()):
            z[i + 1] = gas_frac * y_i

        # Create flash calculator
        # For gamma_phi and explicit: use freshwater flash (salinity=0)
        # For embedded: use salinity in kij correlations directly
        if salinity_method == 'embedded':
            flash = SWMultiComponentFlash(comp_names,
                                           salinity_molal=salinity_molal,
                                           framework=framework)
            result = flash.calc_equilibrium(T_K, P_Pa, z,
                                             salinity_method='embedded')
        else:
            flash = SWMultiComponentFlash(comp_names,
                                           salinity_molal=salinity_molal,
                                           framework=framework)
            result = flash.calc_equilibrium(T_K, P_Pa, z,
                                             salinity_method=salinity_method)

        x_gas = {}
        for i, name in enumerate(comp_names):
            if name != 'H2O':
                if salinity_method == 'explicit' and 'x_aq_brine' in result:
                    x_gas[name] = result['x_aq_brine'][i]
                else:
                    # gamma_phi: x_aq already includes salinity effect
                    # embedded: x_aq already includes salinity effect
                    x_gas[name] = result['x_aq'][i]

        # Water content from kij_NA flash
        y_H2O = result['y_na'][0]  # Water mole fraction from NA flash

    else:
        # Henry's law (binary decomposition) — original method
        x_gas = {}
        for gas, y_i in gas_comp.items():
            vle = SWBinaryVLE(gas, salinity_molal, framework=framework)
            x_i_pure = vle.calc_gas_solubility(T_K, P_Pa)
            x_gas[gas] = x_i_pure * y_i

        # Water content using composition-weighted kij_NA
        kij_na_weighted = sum(y_i * get_kij_na(gas, T_K)
                              for gas, y_i in gas_comp.items())
        dominant_gas = max(gas_comp, key=gas_comp.get)
        vle = SWBinaryVLE(dominant_gas, salinity_molal, framework=framework)
        y_H2O = vle.calc_water_content_with_kij(T_K, P_Pa, kij_na_weighted)

    # Build water content result dict with multiple units
    water_content = {
        'y_H2O': y_H2O,
        'stb_mmscf': y_h2o_to_stb_mmscf(y_H2O),
        'lb_mmscf': y_h2o_to_lb_mmscf(y_H2O),
    }

    return x_gas, water_content


# =============================================================================
# Backward-Compatible Aliases
# =============================================================================
# Alias for existing code that uses BinaryVLE
BinaryVLE = SWBinaryVLE

# Convenience factory for H2-water VLE (matches existing generate_figures.py usage)
class H2WaterVLE(SWBinaryVLE):
    """
    H2-Water VLE calculator. Convenience wrapper for SWBinaryVLE('H2', ...).

    Provides backward compatibility with existing code in generate_figures.py
    and fit_aqueous_bip.py.
    """
    def __init__(self, salinity: float = 0.0):
        """
        Initialize H2-Water VLE calculator.

        Args:
            salinity: Salt concentration in mol/kg water (default: 0 = fresh water)
        """
        super().__init__('H2', salinity)


# Correlation aliases for backward compatibility
def kij_aq_rational(T_K: float) -> float:
    """Rational form kij_AQ for H2 (alias for kij_aq_h2)."""
    return kij_aq_h2(T_K, 0.0)


def kij_aq_linear(T_K: float) -> float:
    """Linear correlation for comparison."""
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for correlation consistency
    return -3.05 + 0.226 * Tr


def kij_na_constant() -> float:
    """Constant non-aqueous phase BIP for H2."""
    return 0.468


def kij_aq_chabab_2023(T_K: float, m: float = 0.0) -> float:
    """
    Chabab et al. 2023 original correlation with integrated salinity effects.

    kij = (D0*(1 + alpha0*m) + D1*Tr*(1 + alpha1*m) + D2*exp(D3*Tr))
    """
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for correlation consistency
    D0, D1, D2, D3 = -2.11917, 0.14888, -13.01835, -0.43946
    alpha0, alpha1 = -0.0226, -0.0045
    return (D0 * (1 + alpha0 * m) + D1 * Tr * (1 + alpha1 * m) + D2 * np.exp(D3 * Tr))


def kij_na_chabab_2023(T_K: float) -> float:
    """Chabab et al. 2023 non-aqueous BIP."""
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for correlation consistency
    return 0.01993 + 0.042834 * Tr


def kij_aq_lopez_lazaro_2019(T_K: float, csw: float = 0.0) -> float:
    """
    Lopez-Lazaro et al. 2019 correlation (Equation 17, Table 6).

    NOTE: Uses exp(-A3*Tr) to avoid non-physical values.
    """
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for correlation consistency
    A0, A1, A2, A3 = -2.513, 0.181, 12.723, 0.499
    a0, a1, b0, b1 = 6.8e-4, 0.038, 0.443, 0.799
    if csw > 0:
        return A0 * (1 + a0 * csw**b0) + A1 * Tr * (1 + a1 * csw**b1) + A2 * np.exp(-A3 * Tr)
    else:
        return A0 + A1 * Tr + A2 * np.exp(-A3 * Tr)


def kij_na_lopez_lazaro_2019(T_K: float) -> float:
    """Lopez-Lazaro et al. 2019 non-aqueous BIP."""
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for correlation consistency
    return 2.500 - 0.179 * Tr


def sechenov_TO(T_K: float, log10_basis: bool = True) -> float:
    """
    Torin-Ollarves & Trusler 2021 Sechenov coefficient (Eq. 22).

    CRITICAL: Original T-O correlation uses NATURAL LOG basis.
    If log10_basis=True (default), convert to log10 basis by dividing by 2.303.
    """
    theta = T_K / 273.15 - 1.0
    d0 = 0.2898
    d1 = -1.4330
    d2 = 3.9584
    d3 = -3.1666
    ks_ln = d0 + d1*theta + d2*theta**2 + d3*theta**3

    if log10_basis:
        return ks_ln / 2.303
    else:
        return ks_ln


# =============================================================================
# Test / Demo
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Soreide-Whitson Unified VLE Engine")
    print("=" * 70)
    print("\nUsing DECOUPLED approach per S&W methodology:")
    print("  - Gas solubility: kij_AQ for both phases")
    print("  - Water content:  kij_NA for both phases")

    # Verify H2 kij_AQ at T=373K
    print("\n" + "=" * 70)
    print("H2 kij_AQ VERIFICATION")
    print("=" * 70)
    T_K = 373.15
    Tr = T_K / BIP_TC_H2  # Use BIP_TC_H2 for kij correlation (33.145 K)
    kij = kij_aq_h2(T_K, 0.0)
    print(f"  T = 373.15 K, Tr = {Tr:.2f} (using BIP_TC_H2 = {BIP_TC_H2} K)")
    print(f"  kij_AQ = ({-14.59} + {Tr:.2f}) / ({2.184} + {0.365}*{Tr:.2f})")
    print(f"         = {-14.59 + Tr:.4f} / {2.184 + 0.365*Tr:.4f}")
    print(f"         = {kij:.4f}")

    # Validation tests
    print("\n" + "=" * 70)
    print("VALIDATION AGAINST EXPERIMENTAL DATA")
    print("=" * 70)

    tests = [
        ('H2', 50, 101.3, 0.001188, 'Wiebe & Gaddy 1934'),
        ('H2', 100, 101.3, 0.001288, 'Wiebe & Gaddy 1934'),
        ('CO2', 50, 100, 0.0163, 'Wiebe & Gaddy 1940'),
        ('CH4', 71, 69, 0.00100, 'Culberson & McKetta'),
        ('N2', 25, 100, 0.00115, 'Wiebe et al. 1933'),
    ]

    print(f"\n{'Gas':<6} {'T(C)':<8} {'P(bar)':<8} {'x_calc':<10} {'x_exp':<10} {'Error%':<10}")
    print("-" * 60)

    for gas, T_C, P_bar, x_exp, source in tests:
        T_K = 273.15 + T_C
        P_Pa = P_bar * 1e5

        vle = SWBinaryVLE(gas, salinity_molal=0.0)
        x_calc = vle.calc_gas_solubility(T_K, P_Pa)
        error = (x_calc - x_exp) / x_exp * 100

        print(f"{gas:<6} {T_C:<8.0f} {P_bar:<8.0f} {x_calc:<10.6f} {x_exp:<10.6f} {error:+8.1f}%")

    # Sechenov verification
    print("\n" + "=" * 70)
    print("SECHENOV COEFFICIENT (S&W Eq. 8)")
    print("=" * 70)
    print(f"\n{'T(C)':<10} {'ks (H2)':<12} {'ks (CH4)':<12} {'ks (N2)':<12}")
    print("-" * 50)
    for T_C in [25, 50, 75, 100, 125, 150]:
        ks_h2 = sw_equation_8_ks(T_C, COMPONENTS['H2'].Tb)
        ks_ch4 = sw_equation_8_ks(T_C, COMPONENTS['CH4'].Tb)
        ks_n2 = sw_equation_8_ks(T_C, COMPONENTS['N2'].Tb)
        print(f"{T_C:<10} {ks_h2:<12.4f} {ks_ch4:<12.4f} {ks_n2:<12.4f}")

    # Multi-component example
    print("\n" + "=" * 70)
    print("MULTI-COMPONENT EXAMPLE")
    print("=" * 70)
    print("\nNatural gas (85% CH4, 10% CO2, 5% N2) in fresh water at 200F, 3000 psia:")

    x_gas, water = calc_gas_brine_equilibrium(
        salinity_wt_pct=0.0,
        temperature_F=200,
        pressure_psia=3000,
        y_CH4=0.85,
        y_CO2=0.10,
        y_N2=0.05
    )

    for gas, x in x_gas.items():
        print(f"  Dissolved {gas}: {x:.6f} mol/mol")
    print(f"  Water in gas: {water['y_H2O']:.6f} mol/mol")
    print(f"               {water['stb_mmscf']:.3f} stb/mmscf")
    print(f"               {water['lb_mmscf']:.1f} lb/mmscf")

    # Gamma-phi flash example
    print("\n" + "=" * 70)
    print("GAMMA-PHI FLASH EXAMPLE (H2 in brine)")
    print("=" * 70)
    print("\n100% H2 in 1.5 molal brine at 50°C, 100 bar:")

    flash = SWMultiComponentFlash(['H2O', 'H2'], salinity_molal=1.5)
    result = flash.calc_equilibrium(323.15, 100e5,
                                     np.array([0.95, 0.05]),
                                     salinity_method='gamma_phi')
    gamma = result.get('gamma', np.ones(2))
    print(f"  γ_H2O = {gamma[0]:.4f}")
    print(f"  γ_H2  = {gamma[1]:.4f}")
    print(f"  x_H2 (brine, gamma-phi) = {result['x_aq'][1]:.6f}")
    print(f"  Converged: {'✓' if result['converged_aq'] else '✗'}")

    # Compare gamma-phi vs explicit Sechenov
    result_expl = flash.calc_equilibrium(323.15, 100e5,
                                          np.array([0.95, 0.05]),
                                          salinity_method='explicit')
    x_expl = result_expl.get('x_aq_brine', result_expl['x_aq'])
    diff = abs(result['x_aq'][1] - x_expl[1]) / result['x_aq'][1] * 100
    print(f"  x_H2 (brine, explicit)  = {x_expl[1]:.6f}")
    print(f"  Difference: {diff:.2f}% (expected <0.3% for dilute H2)")

    # Robust RR solver verification
    print("\n" + "=" * 70)
    print("ROBUST RR SOLVER (Nielsen & Lia 2022)")
    print("=" * 70)
    z_test = np.array([0.90, 0.05, 0.03, 0.02])
    K_test = np.array([0.001, 50.0, 80.0, 120.0])
    V_rr, x_rr, y_rr = solve_rachford_rice(z_test, K_test)
    print(f"\n  z = {z_test}")
    print(f"  K = {K_test}")
    print(f"  V = {V_rr:.6f}")
    print(f"  x = {x_rr}")
    print(f"  y = {y_rr}")
    print(f"  Mass balance: Σ(z-Vy-(1-V)x) = {np.max(np.abs(z_test - V_rr*y_rr - (1-V_rr)*x_rr)):.2e}")

    print("\n" + "=" * 70)
    print("Done.")
