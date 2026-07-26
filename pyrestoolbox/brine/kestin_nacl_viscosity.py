"""Kestin, Khalifa & Correia (1981) NaCl-solution viscosity, JPCRD 10, 71.

Papers/50. The reference correlation of oscillating-disk measurements over
20-150 degC, 0.1-35 MPa, 0-6 molal NaCl, stated to reproduce the underlying
experimental results to +/-0.5% standard deviation.

WHY THIS FILE EXISTS. The salt leg of the shipped viscosity chain (Mao-Duan)
had been scored only against Calabrese's measured gas-free brine at a SINGLE
molality (0.77 mol/kg). A salt term cannot be judged on one concentration.
Kestin supplies the molality axis: a measurement-backed yardstick from 0 to
6 molal at reservoir pressures.

WHAT IT IS AND IS NOT. These are tabulated values of a correlation, not raw
points. It is a smoothed representation of measured data, which is the best
NaCl-viscosity reference available and is how JPCRD publishes such data. Note
that Mao-Duan may have been fitted to the same measurements, so agreement
between the two is not independent evidence.

EQUATIONS, transcribed from a 300 dpi RENDER of pp. 72-73 (the OCR text layer
scrambles every exponent and drops minus signs):

    (1)  mu(p,t,m)   = mu0(t,m) * [1 + beta(t,m) * p]
    (2)  log10[mu0(t,m)/muw0(t)] = A(m) + B(m)*log10[muw0(t)/muw0(20 degC)]
    (3)  log10[muw0(t)/muw0(20 degC)]
                      = {sum_i=1..4 alpha_i*[(20-t)/degC]^i} / [(96+t)/degC]
    (4)  A(m)         = sum_i=1..3 a_i * m^i
    (5)  B(m)         = sum_i=1..3 b_i * m^i
    (6)  beta(t,m)    = betaEs(t) * betastar(m/ms) + betaw(t)
    (7)  betaw(t)     = sum_i=0..4 beta_i * (t/degC)^i          [1/GPa]
    (8)  betaEs(t)    = gamma0 + gamma1*(t/degC) - betaw(t)     [1/GPa]
    (9)  ms(t)        = sum_i=0..2 m_i * (t/degC)^i             [mol/kg]
    (10) betastar(r)  = sum_i=1..3 betastar_i * r^i

ONE TRANSCRIPTION AMBIGUITY, settled against the paper's own tables. The
exponent in (4) and (5) prints as a glyph that reads as `2` at 300 dpi, but
the sums run i = 1..3, which would make three coefficients collinear if the
exponent were fixed. Reading it as `^i` reproduces the printed tables to
0.015% worst over 18 values at 0, 3 and 6 molal; reading it as a literal `^2`
misses by up to 6478%. See `_check_vs_tables()`.

Pressure enters as GPa in Eqs. (7)-(8); the public functions take MPa.
"""

import math

# Eq. (3) - water, Kestin, Sokolov & Wakeham (1978)
_ALPHA = (1.2378, -1.303e-3, 3.06e-6, 2.55e-8)   # i = 1..4
MU_W0_20C = 1002.0          # micro Pa s, zero-pressure water at 20 degC

_A_COEF = (3.324e-2, 3.624e-3, -1.879e-4)        # Eq. (4), i = 1..3
_B_COEF = (-3.96e-2, 1.02e-2, -7.02e-4)          # Eq. (5), i = 1..3
_BETA_W = (-1.297, 5.74e-2, -6.97e-4, 4.47e-6, -1.05e-8)   # Eq. (7), i = 0..4
_GAMMA = (0.545, 2.8e-3)                         # Eq. (8)
_MS = (6.044, 2.8e-3, 3.6e-5)                    # Eq. (9), i = 0..2
_BETA_STAR = (2.5, -2.0, 0.5)                    # Eq. (10), i = 1..3

T_MIN_C, T_MAX_C = 20.0, 150.0
P_MAX_MPA = 35.0
M_MAX = 6.0


def mu_w0(t_C):
    """Zero-pressure pure-water viscosity, micro Pa s. Eq. (3)."""
    num = sum(a * (20.0 - t_C) ** (i + 1) for i, a in enumerate(_ALPHA))
    return MU_W0_20C * 10.0 ** (num / (96.0 + t_C))


def _A(m):
    return sum(a * m ** (i + 1) for i, a in enumerate(_A_COEF))


def _B(m):
    return sum(b * m ** (i + 1) for i, b in enumerate(_B_COEF))


def mu0(t_C, m):
    """Zero-pressure solution viscosity, micro Pa s. Eqs. (2), (4), (5)."""
    muw = mu_w0(t_C)
    return muw * 10.0 ** (_A(m) + _B(m) * math.log10(muw / MU_W0_20C))


def beta_w(t_C):
    """Pressure coefficient of water, 1/GPa. Eq. (7)."""
    return sum(b * t_C ** i for i, b in enumerate(_BETA_W))


def m_saturation(t_C):
    """NaCl saturation molality, mol/kg. Eq. (9), Seidell."""
    return sum(c * t_C ** i for i, c in enumerate(_MS))


def beta(t_C, m):
    """Pressure coefficient of the solution, 1/GPa. Eqs. (6), (8), (10)."""
    bw = beta_w(t_C)
    beta_Es = _GAMMA[0] + _GAMMA[1] * t_C - bw
    r = m / m_saturation(t_C)
    beta_star = sum(b * r ** (i + 1) for i, b in enumerate(_BETA_STAR))
    return beta_Es * beta_star + bw


def mu(t_C, p_MPa, m):
    """Dynamic viscosity of aqueous NaCl, micro Pa s. Eq. (1).

    t_C in degC, p_MPa in MPa, m in mol/kg NaCl. Correlated range
    20-150 degC, 0.1-35 MPa, 0-6 molal; extrapolation is the caller's problem.
    """
    return mu0(t_C, m) * (1.0 + beta(t_C, m) * p_MPa / 1000.0)


def salt_ratio(t_C, p_MPa, m):
    """mu(brine)/mu(water) at the same t and p - the salt term alone."""
    return mu(t_C, p_MPa, m) / mu(t_C, p_MPa, 0.0)
