"""
Pure water properties for the Plyasunov model.

Provides:
    - rho_w(T, P): pure water density in kg/m3
    - kappa_T(T, P): isothermal compressibility in 1/MPa

Backend: IAPWS-IF97 Region 1 (no external dependencies).

Parameters:
    T: temperature in Kelvin
    P: pressure in MPa
"""

from .iapws_if97 import rho_if97, kappa_T_if97

# Constants
MW_WATER = 18.015268   # g/mol, molar mass of water (IAPWS)
TC_WATER = 647.096     # K, critical temperature of water


def rho_w(T, P):
    """
    Pure water density.

    Parameters:
        T: temperature in K
        P: pressure in MPa

    Returns:
        density in kg/m3
    """
    return rho_if97(float(T), float(P))


def kappa_T(T, P):
    """
    Isothermal compressibility of pure water.

    Parameters:
        T: temperature in K
        P: pressure in MPa

    Returns:
        kappa_T in 1/MPa
    """
    return kappa_T_if97(float(T), float(P))
