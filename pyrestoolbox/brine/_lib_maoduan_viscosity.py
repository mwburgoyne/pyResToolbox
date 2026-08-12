"""Mao-Duan (2009) water and NaCl viscosity, the SUPERSEDED route.

Ported from brine_props `salt_viscosity_benchmark.py`. Retained only so
`viscosity_route(route='mao_duan')` reproduces pre-3.7.3 results (to within
0.0094% since the d1/b2 fix below); it is NaCl
only and carries no pressure dependence on the salt ratio. Coefficients are
Mao & Duan (2009) Tables 2-3, verified against the primary 2026-08-08; the
McCain Table 4-14 lineage carried two last-digit defects (d1, b2, <=0.0094%).

Against Kestin's measured NaCl viscosity this ties the default Jones-Dole leg
(0.443% against 0.415% mean) and loses at reservoir temperature (0.754% against
0.408% over 125-150 degC). It was replaced on CAPABILITY: forced to treat KCl as
NaCl it scores 15.9% mean / 51.2% max.
"""
import math

from ..plyasunov.iapws_if97 import rho_if97

_MD_D = [0, 2885317.0, -11072.577, -9.0834095, 0.030925651, -0.0000274071,
         -1928385.1, 5621.6046, 13.82725, -0.047609523, 0.000035545041]
_MD_A = [-0.21319213, 0.0013651589, -0.0000012191756]
_MD_B = [0.069161945, -0.00027292263, 0.00000020852448]
_MD_C = [-0.0025988855, 0.0000077989227]


def mu_water_maoduan(T_K, P_MPa):
    """Mao-Duan pure-water viscosity, mPa s. Density from IF97, as shipped."""
    rho = rho_if97(T_K, P_MPa) / 1000.0
    ln_mu = sum(_MD_D[i] * T_K ** (i - 3) for i in range(1, 6))
    ln_mu += sum(rho * _MD_D[i] * T_K ** (i - 8) for i in range(6, 11))
    return math.exp(ln_mu) * 1e3


def maoduan_ratio(T_K, m):
    """Mao-Duan salt ratio, Eqs. 4.43-4.47. No pressure dependence."""
    A = _MD_A[0] + _MD_A[1] * T_K + _MD_A[2] * T_K * T_K
    B = _MD_B[0] + _MD_B[1] * T_K + _MD_B[2] * T_K * T_K
    C = _MD_C[0] + _MD_C[1] * T_K
    return math.exp(A * m + B * m * m + C * m ** 3)
