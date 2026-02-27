"""
IAPWS-IF97 Region 1: Compressed liquid water properties.

Standalone implementation with no external dependencies beyond math/numpy.
Designed for portability to Dart/Flutter (ResToolbox3).

Provides:
    - rho_if97(T, P): pure water density in kg/m3
    - kappa_T_if97(T, P): isothermal compressibility in 1/MPa

Valid range (Region 1):
    273.15 K <= T <= 623.15 K  (0-350 C)
    P_sat(T) <= P <= 100 MPa

Reference:
    Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997
    for the Thermodynamic Properties of Water and Steam."
    ASME J. Eng. Gas Turbines Power, 122(1), 150-182.

Units: T in K, P in MPa
"""

import numpy as np

# ============================================================================
# Constants
# ============================================================================

R_SPECIFIC = 461.526e-6  # Specific gas constant for water [MPa*m3/(kg*K)]
                         # = 461.526 J/(kg*K) = 0.461526 kJ/(kg*K)
P_STAR = 16.53           # Reference pressure [MPa]
T_STAR = 1386.0          # Reference temperature [K]

# ============================================================================
# Region 1 coefficients (Table 2 of IAPWS-IF97)
# Each row: (I_i, J_i, n_i)
# ============================================================================

_REGION1_IJN = [
    (0,  -2,   0.14632971213167e+00),
    (0,  -1,  -0.84548187389013e+00),
    (0,   0,  -0.37563603672040e+01),
    (0,   1,   0.33855169168385e+01),
    (0,   2,  -0.95791963387872e+00),
    (0,   3,   0.15772038513228e+00),
    (0,   4,  -0.16616417199501e-01),
    (0,   5,   0.81214629983568e-03),
    (1,  -9,   0.28319080123804e-03),
    (1,  -7,  -0.60706301565874e-03),
    (1,  -1,  -0.18990068218419e-01),
    (1,   0,  -0.32529748770505e-01),
    (1,   1,  -0.21841717175414e-01),
    (1,   3,  -0.52838357969930e-04),
    (2,  -3,  -0.47184321073267e-03),
    (2,   0,  -0.30001780793026e-03),
    (2,   1,   0.47661393906987e-04),
    (2,   3,  -0.44141845330846e-05),
    (2,  17,  -0.72694996297594e-15),
    (3,  -4,  -0.31679644845054e-04),
    (3,   0,  -0.28270797985312e-05),
    (3,   6,  -0.85205128120103e-09),
    (4,  -5,  -0.22425281908000e-05),
    (4,  -2,  -0.65171222895601e-06),
    (4,  10,  -0.14341729937924e-12),
    (5,  -8,  -0.40516996860117e-06),
    (8, -11,  -0.12734301741682e-08),
    (8,  -6,  -0.17424871230634e-09),
    (21, -29, -0.68762131295531e-18),
    (23, -31,  0.14478307828521e-19),
    (29, -38,  0.26335781662795e-22),
    (30, -39, -0.11947622640071e-22),
    (31, -40,  0.18228094581404e-23),
    (32, -41, -0.93537087292458e-25),
]


def _gamma_derivatives(pi, tau):
    """
    Compute gamma_pi and gamma_pipi for Region 1.

    gamma = sum( n_i * (7.1 - pi)^I_i * (tau - 1.222)^J_i )

    gamma_pi  = sum( n_i * (-I_i) * (7.1 - pi)^(I_i - 1) * (tau - 1.222)^J_i )
    gamma_pipi = sum( n_i * I_i*(I_i-1) * (7.1 - pi)^(I_i - 2) * (tau - 1.222)^J_i )

    Returns:
        (gamma_pi, gamma_pipi)
    """
    a = 7.1 - pi
    b = tau - 1.222

    gp = 0.0
    gpp = 0.0

    for I, J, n in _REGION1_IJN:
        bJ = b ** J
        if I == 0:
            # gamma_pi term = 0 (I=0 means no pi dependence in derivative)
            # gamma_pipi term = 0
            continue
        elif I == 1:
            # gamma_pi: n * (-1) * (7.1-pi)^0 * bJ = -n * bJ
            gp += -n * bJ
            # gamma_pipi: n * 1*0 * ... = 0
        else:
            aI1 = a ** (I - 1)
            gp += n * (-I) * aI1 * bJ
            gpp += n * I * (I - 1) * (aI1 / a) * bJ  # a^(I-2) = a^(I-1)/a

    return gp, gpp


def rho_if97(T, P):
    """
    Pure water density from IAPWS-IF97 Region 1.

    Parameters:
        T: temperature in K (273.15 - 623.15)
        P: pressure in MPa (up to 100)

    Returns:
        density in kg/m3
    """
    pi = P / P_STAR
    tau = T_STAR / T

    gp, _ = _gamma_derivatives(pi, tau)

    # Specific volume: v = R_SPECIFIC * T * gamma_pi / P_STAR  [m3/kg]
    # Density: rho = P_STAR / (R_SPECIFIC * T * gamma_pi)  [kg/m3]
    return P_STAR / (R_SPECIFIC * T * gp)


def kappa_T_if97(T, P):
    """
    Isothermal compressibility of pure water from IAPWS-IF97 Region 1.

    Parameters:
        T: temperature in K (273.15 - 623.15)
        P: pressure in MPa (up to 100)

    Returns:
        kappa_T in 1/MPa
    """
    pi = P / P_STAR
    tau = T_STAR / T

    gp, gpp = _gamma_derivatives(pi, tau)

    # kappa_T = -gamma_pipi / (gamma_pi * P_STAR)  [1/MPa]
    return -gpp / (gp * P_STAR)


def rho_and_kappa(T, P):
    """
    Compute both density and compressibility in a single call (avoids
    redundant polynomial evaluation).

    Parameters:
        T: temperature in K
        P: pressure in MPa

    Returns:
        (rho [kg/m3], kappa_T [1/MPa])
    """
    pi = P / P_STAR
    tau = T_STAR / T

    gp, gpp = _gamma_derivatives(pi, tau)

    rho = P_STAR / (R_SPECIFIC * T * gp)
    kappa = -gpp / (gp * P_STAR)

    return rho, kappa


# ============================================================================
# Validation
# ============================================================================

if __name__ == "__main__":
    print("=== IAPWS-IF97 Region 1 Validation ===\n")

    # IF97 verification table values (from IAPWS-IF97 document Table 5)
    # T=300K, P=3 MPa:  v = 0.100215168e-2 m3/kg  -> rho = 997.853
    # T=300K, P=80 MPa: v = 0.971180894e-3 m3/kg  -> rho = 1029.672
    # T=500K, P=3 MPa:  v = 0.120241800e-2 m3/kg  -> rho = 831.658
    print("IF97 verification table test points:")
    test_points = [
        (300.0, 3.0,  0.100215168e-2),
        (300.0, 80.0, 0.971180894e-3),
        (500.0, 3.0,  0.120241800e-2),
    ]

    for T, P, v_expected in test_points:
        rho = rho_if97(T, P)
        v = 1.0 / rho
        err = abs(v - v_expected) / v_expected * 100
        print(f"  T={T:.0f}K, P={P:.0f}MPa: v={v:.10e} (expected {v_expected:.10e}) err={err:.6f}%")

    # Compare against iapws library
    print("\n\nComparison with iapws library (IAPWS-95):")
    try:
        from iapws import IAPWS95

        print(f"\n  {'T(C)':>6} {'P(MPa)':>7} {'IF97 rho':>12} {'IAPWS95 rho':>12} {'diff':>10} {'IF97 kT':>12} {'IAPWS95 kT':>12} {'diff':>10}")
        print("  " + "-" * 95)

        conditions = [
            (298.15, 0.1), (298.15, 10.0), (298.15, 50.0),
            (323.15, 10.0), (373.15, 10.0), (373.15, 30.0),
            (423.15, 30.0), (473.15, 30.0), (523.15, 30.0),
            (573.15, 30.0), (373.15, 100.0), (523.15, 100.0),
        ]

        max_rho_err = 0
        max_kt_err = 0

        for T, P in conditions:
            rho_97 = rho_if97(T, P)
            kt_97 = kappa_T_if97(T, P)

            w = IAPWS95(T=T, P=P)
            rho_95 = w.rho
            kt_95 = w.kappa  # 1/MPa

            rho_err = abs(rho_97 - rho_95)
            kt_err = abs(kt_97 - kt_95) / kt_95 * 100

            max_rho_err = max(max_rho_err, rho_err)
            max_kt_err = max(max_kt_err, kt_err)

            T_C = T - 273.15
            print(f"  {T_C:6.1f} {P:7.1f} {rho_97:12.4f} {rho_95:12.4f} {rho_97-rho_95:+10.4f} {kt_97:12.6e} {kt_95:12.6e} {(kt_97-kt_95)/kt_95*100:+9.4f}%")

        print(f"\n  Max density difference: {max_rho_err:.4f} kg/m3")
        print(f"  Max kappa_T relative error: {max_kt_err:.4f}%")

    except ImportError:
        print("  iapws library not available for comparison")
