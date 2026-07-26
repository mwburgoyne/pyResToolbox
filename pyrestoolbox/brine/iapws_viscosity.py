"""
IAPWS-2008 viscosity of ordinary water substance, implemented from scratch.

Reference:
    Huber, M.L. et al. (2009). "New International Formulation for the Viscosity
    of H2O." J. Phys. Chem. Ref. Data 38(2), 101-125.  (Papers/41)

What is implemented: the IAPWS recommendation for industrial use, Huber Eq. (36),

    mu = mu_0(Tbar) * mu_1(Tbar, rhobar) * mu*

which is Eq. (2) with the critical enhancement mu_2 set to 1. Coefficients are
Huber Table 2 (mu_0, 4 terms) and Table 3 (mu_1, 21 terms), both transcribed from
rendered page images rather than the PDF text layer.

What is NOT implemented, and why: the critical enhancement mu_2 (Eqs. 20-33)
contributes more than the 2% uncertainty of the correlation ONLY inside
645.91 K < T < 650.77 K and 245.8 < rho < 405.3 kg/m3 (Huber Eq. 34). That box
lies 200 K above this project's vouched accuracy ceiling of 450 K and entirely
outside IAPWS-IF97 Region 1 (T_max 623.15 K), so no state this code can reach
needs it. It is also unverifiable here: Huber Table 7 requires (d rhobar/d pbar)
from IAPWS-95 at both T and 1.5*T_c, which no dependency-free module on disk
supplies. Shipping unverified code was the worse option.

Also provides Huber Eq. (37) with Table 8 (Patek et al.), the simplified
correlation for liquid water at 0.1 MPa over 253.15-383.15 K.

No external dependencies beyond math, so this ports directly to Dart/Flutter.

Units: T in K, rho in kg/m3, P in MPa, viscosity returned in Pa s.
"""

import math

# ============================================================================
# Reference constants - Huber Eq. (1), p. 104
# ============================================================================

T_STAR = 647.096      # K, critical temperature
RHO_STAR = 322.0      # kg/m3, critical density
P_STAR = 22.064       # MPa, critical pressure
MU_STAR = 1.0e-6      # Pa s

# Range over which the formulation is defined for the liquid, from Huber Sec. 4.
# The project uses it only inside IF97 Region 1, which is narrower.
T_MIN = 253.15
T_MAX = 1173.15

# ============================================================================
# Huber Table 2 (p. 109) - coefficients H_i in Eq. (11) for mu_0(Tbar)
# ============================================================================

_H = [
    1.677_52,
    2.204_62,
    0.636_656_4,
    -0.241_605,
]

# ============================================================================
# Huber Table 3 (p. 111) - coefficients H_ij in Eq. (12) for mu_1(Tbar, rhobar).
# 21 non-zero terms; every H_ij omitted from the table is identically zero.
# Stored as (i, j, H_ij) in the paper's own row order so the table can be
# checked against the page image line by line.
# ============================================================================

_HIJ = [
    (0, 0,  5.200_94e-1),
    (1, 0,  8.508_95e-2),
    (2, 0, -1.083_74),
    (3, 0, -2.895_55e-1),
    (0, 1,  2.225_31e-1),
    (1, 1,  9.991_15e-1),
    (2, 1,  1.887_97),
    (3, 1,  1.266_13),
    (5, 1,  1.205_73e-1),
    (0, 2, -2.813_78e-1),
    (1, 2, -9.068_51e-1),
    (2, 2, -7.724_79e-1),
    (3, 2, -4.898_37e-1),
    (4, 2, -2.570_40e-1),
    (0, 3,  1.619_13e-1),
    (1, 3,  2.573_99e-1),
    (0, 4, -3.253_72e-2),
    (3, 4,  6.984_52e-2),
    (4, 5,  8.721_02e-3),
    (3, 6, -4.356_73e-3),
    (5, 6, -5.932_64e-4),
]

# ============================================================================
# Huber Table 8 (p. 116) - coefficients a_i, b_i in Eq. (37), liquid water at
# 0.1 MPa (Patek et al.). Tbar here is T/(300 K), not T/T_c.
# ============================================================================

_PATEK_AB = [
    (280.68, -1.9),
    (511.45, -7.7),
    (61.131, -19.6),
    (0.459_03, -40.0),
]

PATEK_T_MIN = 253.15
PATEK_T_MAX = 383.15


# ============================================================================
# The correlation
# ============================================================================

def mu0_bar(T_K):
    """
    Dimensionless viscosity in the zero-density limit, Huber Eq. (11):

        mu_0 = 100 * sqrt(Tbar) / sum_{i=0}^{3} H_i / Tbar^i

    Function of temperature only.
    """
    Tbar = T_K / T_STAR
    denom = 0.0
    Tpow = 1.0
    for Hi in _H:
        denom += Hi / Tpow
        Tpow *= Tbar
    return 100.0 * math.sqrt(Tbar) / denom


def mu1_bar(T_K, rho):
    """
    Dimensionless residual (density-dependent) contribution, Huber Eq. (12):

        mu_1 = exp[ rhobar * sum_{i=0}^{5} (1/Tbar - 1)^i
                                 * sum_{j=0}^{6} H_ij (rhobar - 1)^j ]
    """
    Tbar = T_K / T_STAR
    rhobar = rho / RHO_STAR

    dT = 1.0 / Tbar - 1.0
    dR = rhobar - 1.0

    total = 0.0
    for i, j, Hij in _HIJ:
        total += Hij * (dT ** i) * (dR ** j)

    return math.exp(rhobar * total)


def mu_iapws2008(T_K, rho):
    """
    Viscosity of water from temperature and density, Huber Eq. (36)
    (the IAPWS recommendation for industrial use, mu_2 = 1).

    Parameters:
        T_K: temperature in K
        rho: density in kg/m3

    Returns:
        viscosity in Pa s
    """
    return mu0_bar(T_K) * mu1_bar(T_K, rho) * MU_STAR


def mu_water_TP(T_K, P_MPa):
    """
    Viscosity of pure liquid water from temperature and pressure, taking the
    density from IAPWS-IF97 Region 1.

    Huber Sec. 3.6 recommends exactly this pairing: IF97 supplies the density
    when the state is fixed by (T, P). Region 1 bounds (273.15-623.15 K, to
    100 MPa) therefore apply and are enforced by iapws_if97.

    Parameters:
        T_K: temperature in K
        P_MPa: pressure in MPa

    Returns:
        viscosity in Pa s
    """
    from ..plyasunov.iapws_if97 import rho_if97
    return mu_iapws2008(T_K, rho_if97(T_K, P_MPa))


def mu_liquid_0p1MPa(T_K):
    """
    Viscosity of liquid water at 0.1 MPa, Huber Eq. (37) with Table 8.

    A 4-term simplified correlation, stated uncertainty 1% in the stable liquid
    region, valid 253.15-383.15 K and explicitly not to be extrapolated. Used
    only as an independent cross-check on the full formulation at atmospheric
    pressure; the delivered path is mu_water_TP.

    Returns:
        viscosity in Pa s
    """
    if not PATEK_T_MIN <= T_K <= PATEK_T_MAX:
        raise ValueError(
            f"Eq. (37) is valid over {PATEK_T_MIN}-{PATEK_T_MAX} K and Huber "
            f"states it must not be extrapolated; got {T_K} K"
        )
    Ttilde = T_K / 300.0
    return sum(a * Ttilde ** b for a, b in _PATEK_AB) * MU_STAR
