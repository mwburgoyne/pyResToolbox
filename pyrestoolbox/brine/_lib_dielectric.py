"""
Bradley & Pitzer (1979) dielectric constant of water and Debye-Huckel slopes.

Reference:
    Bradley, D.J. and Pitzer, K.S. (1979). "Thermodynamics of electrolytes. 12.
    Dielectric properties of water and Debye-Huckel parameters to 350 degC and
    1 kbar." J. Phys. Chem. 83(12), 1599-1603.  DOI 10.1021/j100475a009
    (Papers/44_BradleyPitzer_1979_dielectric_DebyeHuckel_slopes.pdf)

This supplies A_V, the Debye-Huckel limiting slope for volume, which is the one
input Rogers & Pitzer's NaCl volumetric model (Papers/38, `rogers_pitzer_nacl`)
does not carry itself - they cite this paper for it as their reference [7].
Taking it from here rather than from Archer or IAPWS keeps A_V and the fitted
U parameters from the SAME source, which is the project rule.

MODEL - dielectric constant, Eqs. (1)-(4), Table I. T in K, P in bar.

    D      = D_1000 + C * ln((B + P) / (B + 1000))          (1)
    D_1000 = U1 * exp(U2*T + U3*T^2)                        (2)
    C      = U4 + U5 / (U6 + T)                             (3)
    B      = U7 + U8/T + U9*T                               (4)

The pressure derivative is analytic, and the paper confirms it in the captions
to Figs. 2 and 3:  (dD/dP)_T = C/(B+P)  and  -(d2D/dP2)_T = C/(B+P)^2.

Stated validity: 0-350 degC, to 2000 bar below 70 degC and 5000 bar above.
That covers the whole of IF97 Region 1, so A_V is defined wherever the rest of
the chain is.

DEBYE-HUCKEL SLOPES (p. 1601)

    A_phi = (1/3) * (2*pi*N0*rho_w/1000)^(1/2) * (e^2/(D*k*T))^(3/2)
    A_V   = -2*A_phi*R*T * [ beta_w - 3*(dlnD/dP)_T ]

**THE PRINTED A_V EQUATION HAS A SIGN ERROR.** The paper prints

    A_V = -2*A_phi*R*T * [ 3*(dlnD/dP)_T + beta_w ]      <- as printed, WRONG

verified at 700 dpi, so it is the paper's error and not a scan artefact. It
cannot be right: both beta_w and (dlnD/dP)_T are positive, so the printed form
makes A_V negative, whereas the slope is positive and Rogers & Pitzer tabulate
it as such. Differentiating their own definition of A_phi settles it -

    dln(A_phi)/dP = (1/2)*beta_w - (3/2)*(dlnD/dP)
    A_V = -4*R*T*(dA_phi/dP)_T                      [Rogers & Pitzer Eq. (18)]
        = -2*A_phi*R*T*[ beta_w - 3*(dlnD/dP) ]

- and the correction is confirmed numerically by `verify_against_table_A1()`,
which scores both readings against Rogers & Pitzer's own tabulated A_V column.
Do not "restore" the printed sign.

Units: T in K, P in bar, rho_w in g/cm3, beta_w in 1/bar, A_V in
(cm3/mol)(kg/mol)^1/2. CGS electrostatic units internally.
"""

import math

from ..plyasunov.iapws_if97 import rho_and_kappa

# ============================================================================
# Table I (p. 1600) - constants U1..U9 in Eqs. (2)-(4)
# ============================================================================

U1 = 3.4279e2
U2 = -5.0866e-3
U3 = 9.4690e-7
U4 = -2.0525
U5 = 3.1159e3
U6 = -1.8289e2
U7 = -8.0325e3
U8 = 4.2142e6
U9 = 2.1417

# ============================================================================
# Physical constants. CGS electrostatic, as the A_phi expression requires.
# Values are the CODATA set current when the paper was written; the modern
# values move A_phi by under 1e-5 relative, checked in __main__.
# ============================================================================

N_AVOGADRO = 6.022045e23      # 1/mol
E_CHARGE_ESU = 4.803242e-10   # esu
K_BOLTZMANN = 1.380662e-16    # erg/K
R_CM3_BAR = 83.1440           # cm3 bar mol^-1 K^-1 (Rogers & Pitzer's value)

T_MIN, T_MAX = 273.15, 623.15  # the paper's 0-350 degC


def dielectric_constant(T_K, P_bar):
    """Static dielectric constant of water, Bradley & Pitzer Eq. (1)."""
    D1000 = U1 * math.exp(U2 * T_K + U3 * T_K ** 2)
    C = U4 + U5 / (U6 + T_K)
    B = U7 + U8 / T_K + U9 * T_K
    return D1000 + C * math.log((B + P_bar) / (B + 1000.0))


def dD_dP(T_K, P_bar):
    """(dD/dP)_T in 1/bar. Analytic; the paper states it in the Fig. 2 caption."""
    C = U4 + U5 / (U6 + T_K)
    B = U7 + U8 / T_K + U9 * T_K
    return C / (B + P_bar)


def dlnD_dP(T_K, P_bar):
    """(d ln D/dP)_T in 1/bar."""
    return dD_dP(T_K, P_bar) / dielectric_constant(T_K, P_bar)


def A_phi(T_K, P_bar, rho_w=None):
    """
    Debye-Huckel slope for the osmotic coefficient, (kg/mol)^1/2.

    rho_w in g/cm3; taken from IF97 Region 1 if not supplied.
    """
    if rho_w is None:
        rho_w = rho_and_kappa(T_K, P_bar / 10.0)[0] / 1000.0
    D = dielectric_constant(T_K, P_bar)
    term1 = math.sqrt(2.0 * math.pi * N_AVOGADRO * rho_w / 1000.0)
    term2 = (E_CHARGE_ESU ** 2 / (D * K_BOLTZMANN * T_K)) ** 1.5
    return term1 * term2 / 3.0


def A_V(T_K, P_bar, rho_w=None, beta_w=None, printed_sign=False):
    """
    Debye-Huckel limiting slope for volume, (cm3/mol)(kg/mol)^1/2.

    beta_w is the isothermal compressibility of water in 1/bar; taken from IF97
    if not supplied.

    printed_sign=True reproduces the paper's printed (erroneous) equation, and
    exists only so verify_against_table_A1() can score it. Never use it.
    """
    if rho_w is None or beta_w is None:
        rho, kappa_per_MPa = rho_and_kappa(T_K, P_bar / 10.0)
        if rho_w is None:
            rho_w = rho / 1000.0
        if beta_w is None:
            beta_w = kappa_per_MPa / 10.0        # 1/MPa -> 1/bar

    aphi = A_phi(T_K, P_bar, rho_w)
    dlnD = dlnD_dP(T_K, P_bar)

    if printed_sign:
        return -2.0 * aphi * R_CM3_BAR * T_K * (3.0 * dlnD + beta_w)
    return -2.0 * aphi * R_CM3_BAR * T_K * (beta_w - 3.0 * dlnD)
