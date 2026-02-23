"""
IAPWS-IF97 Region 1: Compressed liquid water properties.

Standalone implementation with no external dependencies beyond numpy.

Provides:
    - rho_if97(T, P): pure water density in kg/m3
    - kappa_T_if97(T, P): isothermal compressibility in 1/MPa
    - rho_and_kappa(T, P): both in a single call

Valid range (Region 1):
    273.15 K <= T <= 623.15 K  (0-350 C)
    P_sat(T) <= P <= 100 MPa

Reference:
    Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997
    for the Thermodynamic Properties of Water and Steam."
    ASME J. Eng. Gas Turbines Power, 122(1), 150-182.

Units: T in K, P in MPa
"""

# Constants
R_SPECIFIC = 461.526e-6  # Specific gas constant for water [MPa*m3/(kg*K)]
P_STAR = 16.53           # Reference pressure [MPa]
T_STAR = 1386.0          # Reference temperature [K]

# Region 1 coefficients (Table 2 of IAPWS-IF97)
# Each row: (I_i, J_i, n_i)
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
            continue
        elif I == 1:
            gp += -n * bJ
        else:
            aI1 = a ** (I - 1)
            gp += n * (-I) * aI1 * bJ
            gpp += n * I * (I - 1) * (aI1 / a) * bJ

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
    return -gpp / (gp * P_STAR)


def rho_and_kappa(T, P):
    """
    Compute both density and compressibility in a single call.

    Returns:
        (rho [kg/m3], kappa_T [1/MPa])
    """
    pi = P / P_STAR
    tau = T_STAR / T
    gp, gpp = _gamma_derivatives(pi, tau)
    rho = P_STAR / (R_SPECIFIC * T * gp)
    kappa = -gpp / (gp * P_STAR)
    return rho, kappa
