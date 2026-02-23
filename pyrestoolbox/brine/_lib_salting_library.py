"""
Library of Sechenov (salting-out) correlations for gases in NaCl brine.

All functions return ks on a LOG10 basis with molality (mol/kg H2O) scale:

    log10(S_water / S_brine) = ks * m_NaCl

Sources implemented:
    1. S&W Eq 8    — Soreide & Whitson 1992, generic (Tb-based), all gases
    2. Dubessy     — Dubessy et al. 2005, extended Sechenov for CO2 and H2S
    3. Akinfiev    — Akinfiev et al. 2016, Pitzer model for H2S (-> effective ks)
    4. Li et al.   — Li, Zhang, Luo & Li 2015, Pitzer for CH4, C2H6, C3H8, nC4H10
    5. Mao & Duan  — Mao & Duan 2006, Pitzer model for N2
    6. Duan & Sun  — Duan & Sun 2003, Pitzer model for CO2

Not implemented (too complex for Sechenov extraction):
    - Springer et al. 2014 (OTC-25295-MS): MSE framework for H2S/CO2.
      Requires speciation + Debye-Huckel + virial + UNIQUAC simultaneously.
      Cannot be reduced to an analytical Sechenov form.

Convention:  ks > 0 means salting-out (lower solubility in brine).
"""

import numpy as np

LN10 = np.log(10.0)  # 2.302585...


# ====================================================================
# 1. Soreide & Whitson 1992, Equation 8  (generic, Tb-based)
# ====================================================================
# Boiling points (K) for S&W gases
TB_K = {
    'H2': 20.3, 'N2': 77.4, 'CO2': 194.7, 'H2S': 213.6,
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
# 2. Dubessy et al. (2005) — Extended Sechenov for CO2 and H2S
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
# 3. Akinfiev, Majer & Shvarov (2016) — Pitzer model for H2S-NaCl
#    Chemical Geology, 424, 1-11. DOI: 10.1016/j.chemgeo.2016.01.006
#
#    Pitzer activity coefficient (Eq. 16):
#      ln(gamma_s) = 2*m_s*lam_ss + 3*m_s^2*tau_sss
#                  + 2*m_e*B_se + 3*m_e^2*C_see + 6*m_s*m_e*C_sse
#
#    For dilute H2S (m_s << m_e), the effective Sechenov is:
#      ks_eff ≈ (2*B_se + 6*m_s*C_sse) / ln(10)   [on log10/molality]
#
#    But the rigorous approach uses the full recommended solubility tables
#    or the complete Pitzer framework.
# ====================================================================

# Binary H2S-H2O self-interaction: Eq. (18)
# lambda_ss = a_lam + b_lam * (100/(T-228)) + c_lam * (T/(T-760))
_AK_A_LAM = -0.19515
_AK_B_LAM = 0.102822
_AK_C_LAM = -0.033726
_AK_TAU_SSS = 0.004900

# Ternary H2S-H2O-NaCl: Eq. (20)
# B_se = b_B * (100/(T-228)) + c_B * (T/(T-760))
_AK_B_B = 0.03568
_AK_C_B = -0.02354
_AK_C_SEE = 0.0       # set to zero by authors
_AK_C_SSE = 0.002558


def _akinfiev_lambda_ss(T_K):
    """H2S-H2S self-interaction parameter (Akinfiev Eq. 18)."""
    T = np.asarray(T_K, dtype=float)
    return _AK_A_LAM + _AK_B_LAM * (100.0 / (T - 228.0)) + _AK_C_LAM * (T / (T - 760.0))


def _akinfiev_B_se(T_K):
    """H2S-NaCl interaction parameter B_se (Akinfiev Eq. 20)."""
    T = np.asarray(T_K, dtype=float)
    return _AK_B_B * (100.0 / (T - 228.0)) + _AK_C_B * (T / (T - 760.0))


def akinfiev_ln_gamma_h2s(T_K, m_h2s, m_NaCl):
    """Natural log of H2S activity coefficient in H2S-H2O-NaCl (Eq. 16).

    Parameters
    ----------
    T_K : float
        Temperature in Kelvin (283-573 K).
    m_h2s : float
        H2S molality (mol/kg H2O).
    m_NaCl : float
        NaCl molality (mol/kg H2O), 0-6.

    Returns
    -------
    ln_gamma : float
        Natural log of activity coefficient on Henry scale.
    """
    lam = _akinfiev_lambda_ss(T_K)
    B = _akinfiev_B_se(T_K)
    return (2.0 * m_h2s * lam + 3.0 * m_h2s**2 * _AK_TAU_SSS
            + 2.0 * m_NaCl * B + 3.0 * m_NaCl**2 * _AK_C_SEE
            + 6.0 * m_h2s * m_NaCl * _AK_C_SSE)


def ks_akinfiev_h2s(T_K, m_NaCl=1.0, m_h2s_approx=0.1):
    """Effective Sechenov coefficient for H2S-NaCl from Akinfiev Pitzer (log10).

    Computes the effective ks from the Pitzer activity coefficient difference
    between saline and fresh solutions, at a specified approximate H2S molality.

    For dilute H2S, ln(gamma_saline/gamma_fresh) ≈ 2*m_e*B_se + 6*m_s*m_e*C_sse
    so ks ≈ (2*B_se + 6*m_s*C_sse) / ln(10).

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


def ks_akinfiev_h2s_from_tables(T_K, m_NaCl):
    """Effective ks for H2S from Akinfiev 2016 recommended solubility tables.

    Uses linearly interpolated Table 3 (pure water) and Table 4 (brine)
    values at 5 MPa to compute ks = log10(m_water/m_brine) / m_NaCl.

    This is the most reliable method — uses the model output directly,
    avoids any approximation in the Pitzer-to-Sechenov conversion.

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin. Supported: 298-523 K.
    m_NaCl : float
        Must be one of: 1, 2, 4, or 6 mol/kg.

    Returns
    -------
    ks : float or array
        Effective Sechenov coefficient (log10 basis).
    """
    # Tables 3 & 4 at P = 5 MPa (mid-range, well within experimental coverage)
    # T_K:                298.15  323.15  373.15  423.15  473.15  523.15
    _T = np.array([298.15, 323.15, 373.15, 423.15, 473.15, 523.15])
    _m0 = np.array([np.nan, 2.026, 1.651, 1.271, 0.956, 0.300])  # pure water, 5 MPa
    # NaCl = 1 m
    _m1 = np.array([np.nan, 1.752, 1.437, 1.110, 0.831, 0.260])
    # NaCl = 2 m
    _m2 = np.array([np.nan, 1.523, 1.259, 0.975, 0.726, 0.226])
    # NaCl = 4 m
    _m4 = np.array([np.nan, 1.168, 0.983, 0.763, 0.561, 0.171])
    # NaCl = 6 m
    _m6 = np.array([np.nan, 0.907, 0.779, 0.606, 0.438, 0.130])

    brine_tables = {1: _m1, 2: _m2, 4: _m4, 6: _m6}
    if m_NaCl not in brine_tables:
        raise ValueError(f"m_NaCl must be 1, 2, 4, or 6. Got {m_NaCl}")

    m_brine = brine_tables[m_NaCl]
    # Compute ks at each table temperature (skip 298.15 where data missing at 5 MPa)
    valid = ~np.isnan(_m0) & ~np.isnan(m_brine)
    ks_pts = np.full_like(_T, np.nan)
    ks_pts[valid] = np.log10(_m0[valid] / m_brine[valid]) / m_NaCl

    # Interpolate to requested temperatures
    T_valid = _T[valid]
    ks_valid = ks_pts[valid]
    T_K = np.asarray(T_K, dtype=float)
    scalar = T_K.ndim == 0
    T_K = np.atleast_1d(T_K)
    result = np.interp(T_K, T_valid, ks_valid)
    return float(result[0]) if scalar else result


# ====================================================================
# 4. Li, Zhang, Luo & Li (2015) — Pitzer model for CH4-C2H6-C3H8-nC4H10
#    Applied Geochemistry. DOI: 10.1016/j.apgeochem.2015.02.006
#
#    Eq. 15 (for NaCl, with lambda_{i-Cl-} = 0):
#      ln(gamma_i) = 2 * m * lambda_{i-Na+}(T,P) + m^2 * zeta_{i-Na+-Cl-}
#
#    Effective Sechenov (log10):
#      ks = (2*lambda + zeta*m) / ln(10)
#
#    Validity: 273-473 K, 1-1000 bar, 0-6 m NaCl
#    NOTE: C3H8 and nC4H10 have very limited brine data (1 atm only).
# ====================================================================

def _li2015_lambda_ch4(T_K, P_bar=100.0):
    """CH4-Na+ Pitzer lambda (Li et al. 2015 Table 7). T in K, P in bar."""
    T = np.asarray(T_K, dtype=float)
    P = float(P_bar)
    return (-5.7066455e-1 + 7.2997588e-4 * T + 1.5176903e2 / T
            + 3.1927112e-5 * P - 1.642651e-5 * P / T)


def _li2015_lambda_c2h6(T_K, P_bar=100.0):
    """C2H6-Na+ Pitzer lambda (Li et al. 2015 Table 7). T in K, P in bar."""
    T = np.asarray(T_K, dtype=float)
    P = float(P_bar)
    return (-2.143686 + 2.598765e-3 * T + 4.6942351e2 / T
            - 4.6849541e-5 * P - 8.4616602e-10 * P**2 * T
            + 1.095219e-6 * P * T)


def _li2015_lambda_c3h8(T_K, P_bar=100.0):
    """C3H8-Na+ Pitzer lambda (Li et al. 2015 Table 7). T-only (no P dep)."""
    T = np.asarray(T_K, dtype=float)
    return 0.513068 - 0.000958 * T


def _li2015_lambda_nc4h10(T_K, P_bar=100.0):
    """nC4H10-Na+ Pitzer lambda (Li et al. 2015 Table 7). T-only (no P dep)."""
    T = np.asarray(T_K, dtype=float)
    return 0.52862384 - 1.0298104e-3 * T


# Zeta (third virial) — all constants, no T or P dependence
_LI2015_ZETA = {
    'CH4':    -2.9990084e-3,
    'C2H6':   -1.0165947e-2,
    'C3H8':   -0.007485,
    'NC4H10':  0.0206946,    # NOTE: positive (unlike the others)
}

_LI2015_LAMBDA = {
    'CH4':    _li2015_lambda_ch4,
    'C2H6':   _li2015_lambda_c2h6,
    'C3H8':   _li2015_lambda_c3h8,
    'NC4H10': _li2015_lambda_nc4h10,
}


def ks_li2015(T_K, gas, m_NaCl=0.0, P_bar=100.0):
    """Effective Sechenov for light HCs from Li et al. 2015 Pitzer (log10).

    ln(gamma) = 2*m*lambda(T,P) + m^2*zeta
    ks = (2*lambda + zeta*m) / ln(10)

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin (273-473 K).
    gas : str
        'CH4', 'C2H6', 'C3H8', or 'nC4H10'.
    m_NaCl : float
        NaCl molality. Affects effective ks through zeta term.
    P_bar : float
        Pressure in bar (affects CH4 and C2H6 lambda). Default 100.

    Returns
    -------
    ks : float or array
        Effective Sechenov coefficient (log10 basis, molality scale).
    """
    g = gas.upper().replace('N', 'N')  # normalize
    if g == 'NC4H10' or g == 'NC4':
        g = 'NC4H10'
    elif g not in _LI2015_LAMBDA:
        raise ValueError(f"Li 2015 covers CH4, C2H6, C3H8, nC4H10. Got {gas}")

    lam = _LI2015_LAMBDA[g](T_K, P_bar)
    zeta = _LI2015_ZETA[g]
    return (2.0 * lam + zeta * m_NaCl) / LN10


# ====================================================================
# 5. Mao & Duan (2006) — Pitzer model for N2-NaCl
#    Fluid Phase Equilibria, 248, 103-114.
#    DOI: 10.1016/j.fluid.2006.07.020
#
#    Eq. 8 (for NaCl, with lambda_{N2-Cl-} = 0):
#      ln(gamma_N2) = 2 * m * lambda_{N2-Na+}(T,P) + m^2 * xi_{N2-Na+-Cl-}
#
#    Effective Sechenov (log10):
#      ks = (2*lambda + xi*m) / ln(10)
#
#    Validity: 273-400 K, 1-600 bar, 0-6 m NaCl
# ====================================================================

def _mao2006_lambda_n2(T_K, P_bar=100.0):
    """N2-Na+ Pitzer lambda (Mao & Duan 2006 Table 3). T in K, P in bar."""
    T = np.asarray(T_K, dtype=float)
    P = float(P_bar)
    return (-2.4434074 + 0.0036351795 * T + 447.47364 / T
            - 0.000013711527 * P + 0.0000071037217 * P**2 / T)


_MAO2006_XI_N2 = -0.0058071053  # constant, no T or P dependence


def ks_mao2006_n2(T_K, m_NaCl=0.0, P_bar=100.0):
    """Effective Sechenov for N2 from Mao & Duan 2006 Pitzer model (log10).

    ln(gamma_N2) = 2*m*lambda(T,P) + m^2*xi
    ks = (2*lambda + xi*m) / ln(10)

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin (273-400 K).
    m_NaCl : float
        NaCl molality. Affects effective ks through xi term.
    P_bar : float
        Pressure in bar (1-600 bar). Default 100.

    Returns
    -------
    ks : float or array
        Effective Sechenov coefficient (log10 basis, molality scale).
    """
    lam = _mao2006_lambda_n2(T_K, P_bar)
    return (2.0 * lam + _MAO2006_XI_N2 * m_NaCl) / LN10


# ====================================================================
# 6. Duan & Sun (2003) — Pitzer model for CO2-NaCl
#    Chemical Geology, 193, 257-271.
#    DOI: 10.1016/S0009-2541(02)00263-2
#
#    Same Pitzer framework as Mao 2006 (N2) — Pitzer et al. (1984) form:
#      Par(T,P) = c1 + c2*T + c3/T + c4*T^2 + c5/(630-T)
#               + c6*P + c7*P*ln(T) + c8*P/T + c9*P/(630-T)
#               + c10*P^2/(630-T)^2 + c11*T*ln(P)
#
#    lambda_{CO2-Na+} uses c1,c2,c3,c8,c9 (5 coeffs)
#    zeta_{CO2-Na+-Cl-} uses c1,c2,c8,c9 (4 coeffs)
#    lambda_{CO2-Cl-} = 0 (set to zero)
#
#    ln(gamma_CO2) = 2*m*lambda + m^2*zeta
#    Effective Sechenov: ks = (2*lambda + zeta*m) / ln(10)
#
#    Validity: 273-533 K, 0-2000 bar, 0-4.3 m NaCl
# ====================================================================

def _duan2003_lambda_co2(T_K, P_bar=100.0):
    """CO2-Na+ Pitzer lambda (Duan & Sun 2003 Table 2). T in K, P in bar."""
    T = np.asarray(T_K, dtype=float)
    P = float(P_bar)
    return (-0.411370585 + 6.07632013e-4 * T + 97.5347708 / T
            - 0.0237622469 * P / T + 0.0170656236 * P / (630.0 - T))


def _duan2003_zeta_co2(T_K, P_bar=100.0):
    """CO2-Na+-Cl- Pitzer zeta (Duan & Sun 2003 Table 2). T in K, P in bar."""
    T = np.asarray(T_K, dtype=float)
    P = float(P_bar)
    return (3.36389723e-4 - 1.98298980e-5 * T
            + 2.12220830e-3 * P / T - 5.24873303e-3 * P / (630.0 - T))


def ks_duan2003_co2(T_K, m_NaCl=0.0, P_bar=100.0):
    """Effective Sechenov for CO2 from Duan & Sun 2003 Pitzer model (log10).

    ln(gamma_CO2) = 2*m*lambda(T,P) + m^2*zeta(T,P)
    ks = (2*lambda + zeta*m) / ln(10)

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin (273-533 K).
    m_NaCl : float
        NaCl molality (0-4.3 m). Affects effective ks through zeta term.
    P_bar : float
        Pressure in bar (0-2000 bar). Default 100.

    Returns
    -------
    ks : float or array
        Effective Sechenov coefficient (log10 basis, molality scale).
    """
    lam = _duan2003_lambda_co2(T_K, P_bar)
    zeta = _duan2003_zeta_co2(T_K, P_bar)
    return (2.0 * lam + zeta * m_NaCl) / LN10


# ====================================================================
# Convenience: unified interface
# ====================================================================

def ks_library(T_K, gas, source='sw_eq8', m_NaCl=0.0, **kwargs):
    """Unified Sechenov coefficient lookup (log10 basis, molality scale).

    Parameters
    ----------
    T_K : float or array
        Temperature in Kelvin.
    gas : str
        Gas name: 'H2', 'N2', 'CH4', 'C2H6', 'C3H8', 'CO2', 'H2S', etc.
    source : str
        'sw_eq8'        — Soreide & Whitson 1992 Eq 8 (all gases)
        'dubessy'       — Dubessy et al. 2005 (CO2, H2S only)
        'akinfiev'      — Akinfiev et al. 2016 Pitzer (H2S only)
        'akinfiev_table'— Akinfiev et al. 2016 from tables (H2S only)
        'li2015'        — Li et al. 2015 Pitzer (CH4, C2H6, C3H8, nC4H10)
        'mao2006'       — Mao & Duan 2006 Pitzer (N2 only)
        'duan2003'      — Duan & Sun 2003 Pitzer (CO2 only)
    m_NaCl : float
        NaCl molality (needed for non-linear Sechenov models).

    Returns
    -------
    ks : float or array
    """
    gas = gas.upper()
    if gas == 'NC4H10':
        pass  # already normalized
    elif gas == 'NC4':
        gas = 'NC4H10'
    src = source.lower()

    if src == 'sw_eq8':
        return ks_sw_eq8(T_K, gas)
    elif src == 'dubessy':
        if gas == 'CO2':
            return ks_dubessy_co2(T_K, m_NaCl)
        elif gas == 'H2S':
            return ks_dubessy_h2s(T_K, m_NaCl)
        else:
            raise ValueError(f"Dubessy only covers CO2 and H2S, not {gas}")
    elif src == 'akinfiev':
        if gas == 'H2S':
            return ks_akinfiev_h2s(T_K, m_NaCl, **kwargs)
        else:
            raise ValueError(f"Akinfiev 2016 only covers H2S, not {gas}")
    elif src == 'akinfiev_table':
        if gas == 'H2S':
            return ks_akinfiev_h2s_from_tables(T_K, m_NaCl)
        else:
            raise ValueError(f"Akinfiev tables only for H2S, not {gas}")
    elif src == 'li2015':
        return ks_li2015(T_K, gas, m_NaCl, **kwargs)
    elif src == 'mao2006':
        if gas == 'N2':
            return ks_mao2006_n2(T_K, m_NaCl, **kwargs)
        else:
            raise ValueError(f"Mao & Duan 2006 only covers N2, not {gas}")
    elif src == 'duan2003':
        if gas == 'CO2':
            return ks_duan2003_co2(T_K, m_NaCl, **kwargs)
        else:
            raise ValueError(f"Duan & Sun 2003 only covers CO2, not {gas}")
    else:
        raise ValueError(f"Unknown source: {source}")


# ====================================================================
# Main: comparison table
# ====================================================================
if __name__ == '__main__':
    temps_C = [25, 50, 75, 100, 125, 150, 200, 250]
    temps_K = [T + 273.15 for T in temps_C]

    print("=" * 80)
    print("Sechenov Coefficient Library — All values on log10 / molality basis")
    print("  ks > 0 means salting-out;  log10(S_water/S_brine) = ks * m_NaCl")
    print("=" * 80)

    # --- CO2 comparison ---
    print("\n--- CO2: S&W vs Dubessy vs Duan 2003 ---")
    print(f"{'T(C)':>6}  {'S&W Eq8':>8}  {'Dub m=0':>8}  {'Dub m=1':>8}  "
          f"{'Dub m=4':>8}  {'Duan m=0':>9}  {'Duan m=1':>9}  {'Duan m=4':>9}")
    print("-" * 82)
    for Tc, Tk in zip(temps_C, temps_K):
        sw = ks_sw_eq8(Tk, 'CO2')
        d0 = ks_dubessy_co2(Tk, 0)
        d1 = ks_dubessy_co2(Tk, 1)
        d4 = ks_dubessy_co2(Tk, 4)
        if Tk <= 533:
            dn0 = ks_duan2003_co2(Tk, 0, 100)
            dn1 = ks_duan2003_co2(Tk, 1, 100)
            dn4 = ks_duan2003_co2(Tk, 4, 100)
            print(f"{Tc:>6}  {sw:8.4f}  {d0:8.4f}  {d1:8.4f}  {d4:8.4f}  "
                  f"{dn0:9.4f}  {dn1:9.4f}  {dn4:9.4f}")
        else:
            print(f"{Tc:>6}  {sw:8.4f}  {d0:8.4f}  {d1:8.4f}  {d4:8.4f}  "
                  f"{'n/a':>9}  {'n/a':>9}  {'n/a':>9}")

    # --- H2S comparison ---
    print("\n--- H2S ---")
    print(f"{'T(C)':>6}  {'S&W Eq8':>8}  {'Dub m=0':>8}  {'Dub m=1':>8}  "
          f"{'Akin m=1':>9}  {'Akin m=4':>9}  {'AkTbl m=1':>10}  {'AkTbl m=4':>10}")
    print("-" * 90)
    for Tc, Tk in zip(temps_C, temps_K):
        sw = ks_sw_eq8(Tk, 'H2S')
        d0 = ks_dubessy_h2s(Tk, 0)
        d1 = ks_dubessy_h2s(Tk, 1)
        ak1 = ks_akinfiev_h2s(Tk, m_NaCl=1.0, m_h2s_approx=0.1)
        ak4 = ks_akinfiev_h2s(Tk, m_NaCl=4.0, m_h2s_approx=0.1)
        # Table-based (only for T >= 323.15 K)
        if Tk >= 323.15:
            akt1 = ks_akinfiev_h2s_from_tables(Tk, 1)
            akt4 = ks_akinfiev_h2s_from_tables(Tk, 4)
            print(f"{Tc:>6}  {sw:8.4f}  {d0:8.4f}  {d1:8.4f}  "
                  f"{ak1:9.4f}  {ak4:9.4f}  {akt1:10.4f}  {akt4:10.4f}")
        else:
            print(f"{Tc:>6}  {sw:8.4f}  {d0:8.4f}  {d1:8.4f}  "
                  f"{ak1:9.4f}  {ak4:9.4f}  {'n/a':>10}  {'n/a':>10}")

    # --- Light HC comparison: S&W vs Li 2015 ---
    print("\n--- Light Hydrocarbons: S&W Eq 8 vs Li et al. 2015 Pitzer ---")
    print("    Li 2015 at P=100 bar, m=0 (infinite dilution) and m=2")
    hc_gases = ['CH4', 'C2H6', 'C3H8', 'nC4H10']
    print(f"{'T(C)':>6}", end="")
    for g in hc_gases:
        print(f"  {'SW '+g:>11}  {'Li m=0':>7}  {'Li m=2':>7}", end="")
    print()
    print("-" * 118)
    for Tc, Tk in zip(temps_C, temps_K):
        print(f"{Tc:>6}", end="")
        for g in hc_gases:
            sw = ks_sw_eq8(Tk, g)
            li0 = ks_li2015(Tk, g, m_NaCl=0.0, P_bar=100.0)
            li2 = ks_li2015(Tk, g, m_NaCl=2.0, P_bar=100.0)
            print(f"  {sw:11.4f}  {li0:7.4f}  {li2:7.4f}", end="")
        print()

    # --- All S&W gases ---
    print("\n--- S&W Eq 8 for all gases (m-independent) ---")
    gases = ['H2', 'N2', 'CH4', 'C2H6', 'C3H8', 'nC4H10', 'CO2', 'H2S']
    header = f"{'T(C)':>6}" + "".join(f"  {g:>8}" for g in gases)
    print(header)
    print("-" * (8 + 10 * len(gases)))
    for Tc, Tk in zip(temps_C, temps_K):
        row = f"{Tc:>6}"
        for g in gases:
            row += f"  {ks_sw_eq8(Tk, g):8.4f}"
        print(row)

    # --- N2 comparison: S&W vs Mao & Duan 2006 ---
    print("\n--- N2: S&W Eq 8 vs Mao & Duan 2006 Pitzer ---")
    print("    Mao 2006 at P=100 bar, various m_NaCl")
    print(f"{'T(C)':>6}  {'S&W Eq8':>8}  {'Mao m=0':>8}  {'Mao m=1':>8}  "
          f"{'Mao m=2':>8}  {'Mao m=4':>8}")
    print("-" * 56)
    for Tc, Tk in zip(temps_C, temps_K):
        if Tk > 400:
            # Mao 2006 valid to 400 K only
            sw = ks_sw_eq8(Tk, 'N2')
            print(f"{Tc:>6}  {sw:8.4f}  {'n/a':>8}  {'n/a':>8}  {'n/a':>8}  {'n/a':>8}")
        else:
            sw = ks_sw_eq8(Tk, 'N2')
            m0 = ks_mao2006_n2(Tk, m_NaCl=0.0, P_bar=100.0)
            m1 = ks_mao2006_n2(Tk, m_NaCl=1.0, P_bar=100.0)
            m2 = ks_mao2006_n2(Tk, m_NaCl=2.0, P_bar=100.0)
            m4 = ks_mao2006_n2(Tk, m_NaCl=4.0, P_bar=100.0)
            print(f"{Tc:>6}  {sw:8.4f}  {m0:8.4f}  {m1:8.4f}  {m2:8.4f}  {m4:8.4f}")

    # --- N2: Pressure sensitivity (Mao 2006) ---
    print("\n--- N2: Mao 2006 ks at different pressures (m=0) ---")
    print(f"{'T(C)':>6}  {'1 bar':>8}  {'50 bar':>8}  {'100 bar':>8}  "
          f"{'200 bar':>8}  {'500 bar':>8}  {'S&W':>8}")
    print("-" * 62)
    for Tc, Tk in zip(temps_C, temps_K):
        if Tk > 400:
            sw = ks_sw_eq8(Tk, 'N2')
            print(f"{Tc:>6}  {'n/a':>8}  {'n/a':>8}  {'n/a':>8}  {'n/a':>8}  {'n/a':>8}  {sw:8.4f}")
        else:
            vals = [ks_mao2006_n2(Tk, 0, P) for P in [1, 50, 100, 200, 500]]
            sw = ks_sw_eq8(Tk, 'N2')
            print(f"{Tc:>6}" + "".join(f"  {v:8.4f}" for v in vals) + f"  {sw:8.4f}")

    # --- CH4: Pressure sensitivity (Li 2015) ---
    print("\n--- CH4: Li 2015 ks at different pressures (m=0) ---")
    print(f"{'T(C)':>6}  {'1 bar':>8}  {'50 bar':>8}  {'100 bar':>8}  "
          f"{'200 bar':>8}  {'500 bar':>8}  {'S&W':>8}")
    print("-" * 62)
    for Tc, Tk in zip(temps_C, temps_K):
        vals = [ks_li2015(Tk, 'CH4', 0, P) for P in [1, 50, 100, 200, 500]]
        sw = ks_sw_eq8(Tk, 'CH4')
        print(f"{Tc:>6}" + "".join(f"  {v:8.4f}" for v in vals) + f"  {sw:8.4f}")

    # --- H2S model comparison summary ---
    print("\n--- H2S: Akinfiev table ks at multiple pressures (log10) ---")
    print("  Showing that ks is largely P-independent (as expected)")
    # Tables at 1 MPa and 10 MPa for m=2
    T_tab = [323.15, 373.15, 423.15, 473.15, 523.15]
    m0_1 = [0.596, 0.297, 0.136, np.nan, np.nan]     # pure water, 1 MPa
    m2_1 = [0.467, 0.241, 0.111, np.nan, np.nan]      # 2m NaCl, 1 MPa
    m0_5 = [2.026, 1.651, 1.271, 0.956, 0.300]        # pure water, 5 MPa
    m2_5 = [1.523, 1.259, 0.975, 0.726, 0.226]        # 2m NaCl, 5 MPa
    m0_10 = [2.100, 2.703, 2.752, 2.622, 2.045]       # pure water, 10 MPa
    m2_10 = [1.577, 1.980, 1.969, 1.834, 1.411]       # 2m NaCl, 10 MPa
    m0_20 = [2.237, 3.065, 4.407, 5.395, 5.458]       # pure water, 20 MPa
    m2_20 = [1.676, 2.224, 3.039, 3.659, 3.677]       # 2m NaCl, 20 MPa

    print(f"{'T(C)':>6}  {'1 MPa':>8}  {'5 MPa':>8}  {'10 MPa':>8}  {'20 MPa':>8}")
    print("-" * 46)
    for i, Tk in enumerate(T_tab):
        Tc = Tk - 273.15
        vals = []
        for mw, mb in [(m0_1[i], m2_1[i]), (m0_5[i], m2_5[i]),
                        (m0_10[i], m2_10[i]), (m0_20[i], m2_20[i])]:
            if np.isnan(mw) or np.isnan(mb) or mb <= 0:
                vals.append("   n/a  ")
            else:
                ks = np.log10(mw / mb) / 2.0
                vals.append(f"{ks:8.4f}")
        print(f"{Tc:>6.0f}  {'  '.join(vals)}")
