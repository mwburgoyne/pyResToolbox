"""
Library of Sechenov (salting-out) correlations for gases in NaCl brine.

All functions return ks on a LOG10 basis with molality (mol/kg H2O) scale:

    log10(S_water / S_brine) = ks * m_NaCl

Models retained (only those routed by the VLE engine, plus reference forms):
    1. S&W Eq 8  - Soreide & Whitson 1992, generic (Tb-based), all gases
    2. Dubessy   - Dubessy et al. 2005, extended Sechenov for CO2 and H2S
    3. Akinfiev  - Akinfiev et al. 2016, Pitzer model for H2S (-> effective ks)

The VLE engine (get_sechenov_ks, proposed framework) routes CO2 to
ks_dubessy_co2 and H2S to ks_akinfiev_h2s; all other gases use S&W Eq 8.

Convention:  ks > 0 means salting-out (lower solubility in brine).
"""

import numpy as np

LN10 = np.log(10.0)  # 2.302585...


# ====================================================================
# 1. Soreide & Whitson 1992, Equation 8  (generic, Tb-based)
# ====================================================================
# Boiling points (K) for S&W gases
TB_K = {
    'H2': 20.3, 'N2': 77.4, 'CO2': 194.7, 'H2S': 212.8,  # H2S NBP 212.87 K (NIST)
    'CH4': 111.6, 'C2H6': 184.6, 'C3H8': 231.1,
    'iC4H10': 261.4, 'nC4H10': 272.7, 'iC5H12': 301.0,
    'nC5H12': 309.2, 'nC6H14': 341.9, 'nC7H16': 371.6,
    'nC8H18': 398.8, 'nC10H22': 447.3,
}


def ks_sw_eq8(T_K, gas_or_Tb):
    """Soreide-Whitson 1992 Equation 8 Sechenov coefficient.

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin.
    gas_or_Tb : str or float
        Gas name (e.g. 'H2', 'CH4') or boiling point in K.

    Returns
    -------
    ks : float or array
        Sechenov coefficient (log10 basis, molality scale).
    """
    if isinstance(gas_or_Tb, str):
        Tb = TB_K[gas_or_Tb]
    else:
        Tb = float(gas_or_Tb)
    T_F = (np.asarray(T_K, dtype=float) - 273.15) * 9.0 / 5.0 + 32.0
    return (0.13163 + 4.45e-4 * Tb - 7.692e-4 * T_F
            + 2.6614e-6 * T_F**2 - 2.612e-9 * T_F**3)


# ====================================================================
# 2. Dubessy et al. (2005) - Extended Sechenov for CO2 and H2S
#    Oil & Gas Sci. Technol. - Rev. IFP, Vol. 60, No. 2, pp. 339-355
#    Table 9 coefficients.
#
#    log10(K_brine/K_water) = m*b1(T) + m^2*b2(T) + m^3*b3(T)
#    b_i(T) = sum(B_ij * T^j, j=0..4)
#
#    Effective ks = b1(T) + m*b2(T) + m^2*b3(T)
# ====================================================================

# CO2: column scaling 1e0, 1e-2, 1e-4, 1e-8, 1e-11
_DUB_CO2_B1 = np.array([3.114712456, -2.7655585e-2, 0.9176713976e-4,
                         -12.78795941e-8, 6.2704268351e-11])
_DUB_CO2_B2 = np.array([-2.05637458, 2.081980200e-2, -0.765857702e-4,
                         12.011325315e-8, -6.790343083e-11])
_DUB_CO2_B3 = np.array([0.253424331, -0.26047432e-2, 0.0972580216e-4,
                         -1.551654794e-8, 0.8948557284e-11])

# H2S: column scaling 1e0, 1e-2, 1e-4, 1e-7, 1e-10
_DUB_H2S_B1 = np.array([-12.4617636, 12.69373100e-2, -4.791540697e-4,
                         7.9817223650e-7, -4.931093145e-10])
_DUB_H2S_B2 = np.array([5.327383011, -5.82779828e-2, 2.3650333285e-4,
                         -4.207913036e-7, 2.7628521914e-10])
_DUB_H2S_B3 = np.array([-0.75715275, 0.831927411e-2, -0.338668040e-4,
                         0.6037602785e-7, -0.397049836e-10])


def _poly4(c, T):
    """Evaluate 4th-order polynomial c[0] + c[1]*T + ... + c[4]*T^4."""
    T = np.asarray(T, dtype=float)
    return c[0] + c[1]*T + c[2]*T**2 + c[3]*T**3 + c[4]*T**4


def ks_dubessy_co2(T_K, m_NaCl=0.0):
    """Dubessy et al. 2005 effective Sechenov for CO2-NaCl (log10).

    Valid: T <= 543 K (270 C), P <= 300 bar, m <= 6. R^2 = 0.83.
    """
    return (_poly4(_DUB_CO2_B1, T_K) + m_NaCl * _poly4(_DUB_CO2_B2, T_K)
            + m_NaCl**2 * _poly4(_DUB_CO2_B3, T_K))


def ks_dubessy_h2s(T_K, m_NaCl=0.0):
    """Dubessy et al. 2005 effective Sechenov for H2S-NaCl (log10).

    Valid: T <= 523 K (250 C), P <= 150 bar, m <= 6. R^2 = 0.74.
    """
    return (_poly4(_DUB_H2S_B1, T_K) + m_NaCl * _poly4(_DUB_H2S_B2, T_K)
            + m_NaCl**2 * _poly4(_DUB_H2S_B3, T_K))


# ====================================================================
# 3. Akinfiev, Majer & Shvarov (2016) - Pitzer model for H2S-NaCl
#    Chemical Geology, 424, 1-11. DOI: 10.1016/j.chemgeo.2016.01.006
#
#    For dilute H2S (m_s << m_e), the effective Sechenov is:
#      ks_eff ~= (2*B_se + 6*m_s*C_sse) / ln(10)   [on log10/molality]
# ====================================================================

# Ternary H2S-H2O-NaCl: Eq. (20)
# B_se = b_B * (100/(T-228)) + c_B * (T/(T-760))
_AK_B_B = 0.03568
_AK_C_B = -0.02354
_AK_C_SEE = 0.0       # set to zero by authors
_AK_C_SSE = 0.002558


def _akinfiev_B_se(T_K):
    """H2S-NaCl interaction parameter B_se (Akinfiev Eq. 20)."""
    T = np.asarray(T_K, dtype=float)
    return _AK_B_B * (100.0 / (T - 228.0)) + _AK_C_B * (T / (T - 760.0))


def ks_akinfiev_h2s(T_K, m_NaCl=1.0, m_h2s_approx=0.1):
    """Effective Sechenov coefficient for H2S-NaCl from Akinfiev Pitzer (log10).

    Computes the effective ks from the Pitzer activity coefficient difference
    between saline and fresh solutions, at a specified approximate H2S molality.

    For dilute H2S, ln(gamma_saline/gamma_fresh) ~= 2*m_e*B_se + 6*m_s*m_e*C_sse
    so ks ~= (2*B_se + 6*m_s*C_sse) / ln(10).

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin (283-573 K).
    m_NaCl : float
        NaCl molality for computing the effective ks. Default 1 m.
    m_h2s_approx : float
        Approximate H2S molality (affects C_sse contribution). Default 0.1 m.

    Returns
    -------
    ks : float or array
        Effective Sechenov coefficient (log10 basis, molality scale).
    """
    B = _akinfiev_B_se(T_K)
    # ln(gamma_brine) - ln(gamma_fresh) at given m_h2s:
    #   = 2*m_NaCl*B_se + 3*m_NaCl^2*C_see + 6*m_h2s*m_NaCl*C_sse
    # Dividing by m_NaCl gives the ln-based Sechenov:
    ks_ln = 2.0 * B + 3.0 * m_NaCl * _AK_C_SEE + 6.0 * m_h2s_approx * _AK_C_SSE
    return ks_ln / LN10
