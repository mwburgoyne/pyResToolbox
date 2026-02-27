"""
Pure water properties for the Plyasunov model.

Provides:
    - rho_w(T, P): pure water density in kg/m3
    - kappa_T(T, P): isothermal compressibility in 1/MPa
    - V1_star(T, P): molar volume of water in cm3/mol

Default backend: IAPWS-IF97 Region 1 (from-scratch, no external dependencies).
Optional: IAPWS-95 via iapws library (set USE_IAPWS95=True).

Parameters:
    T: temperature in Kelvin
    P: pressure in MPa
"""

import numpy as np
from iapws_if97 import rho_if97, kappa_T_if97, rho_and_kappa

# Set to True to use the full IAPWS-95 formulation via iapws library
USE_IAPWS95 = False

# Constants
MW_WATER = 18.015268  # g/mol, molar mass of water (IAPWS)
TC_WATER = 647.096     # K, critical temperature of water
PC_WATER = 22.064      # MPa, critical pressure of water
RHOC_WATER = 322.0     # kg/m3, critical density of water
V1C_STAR = 55.95       # cm3/mol, critical molar volume of water


def rho_w(T, P):
    """
    Pure water density.

    Parameters:
        T: temperature in K (scalar or array)
        P: pressure in MPa (scalar or array)

    Returns:
        density in kg/m3
    """
    if USE_IAPWS95:
        return _rho_w_iapws95(T, P)

    T = np.atleast_1d(np.asarray(T, dtype=float))
    P = np.atleast_1d(np.asarray(P, dtype=float))
    T, P = np.broadcast_arrays(T, P)

    result = np.zeros_like(T)
    for i in np.ndindex(T.shape):
        result[i] = rho_if97(float(T[i]), float(P[i]))

    return float(result) if result.ndim == 0 else result.squeeze()


def kappa_T(T, P):
    """
    Isothermal compressibility of pure water.

    Parameters:
        T: temperature in K
        P: pressure in MPa

    Returns:
        kappa_T in 1/MPa
    """
    if USE_IAPWS95:
        return _kappa_T_iapws95(T, P)

    T = np.atleast_1d(np.asarray(T, dtype=float))
    P = np.atleast_1d(np.asarray(P, dtype=float))
    T, P = np.broadcast_arrays(T, P)

    result = np.zeros_like(T)
    for i in np.ndindex(T.shape):
        result[i] = kappa_T_if97(float(T[i]), float(P[i]))

    return float(result) if result.ndim == 0 else result.squeeze()


def V1_star(T, P):
    """
    Molar volume of pure water in cm3/mol.

    Parameters:
        T: temperature in K
        P: pressure in MPa

    Returns:
        V1* in cm3/mol
    """
    rho = rho_w(T, P)
    return MW_WATER * 1000.0 / rho


# ============================================================================
# IAPWS-95 backend (requires iapws library)
# ============================================================================

def _rho_w_iapws95(T, P):
    from iapws import IAPWS95
    T = np.atleast_1d(np.asarray(T, dtype=float))
    P = np.atleast_1d(np.asarray(P, dtype=float))
    T, P = np.broadcast_arrays(T, P)
    result = np.zeros_like(T)
    for i in np.ndindex(T.shape):
        w = IAPWS95(T=float(T[i]), P=float(P[i]))
        result[i] = w.rho
    return float(result) if result.ndim == 0 else result.squeeze()


def _kappa_T_iapws95(T, P):
    from iapws import IAPWS95
    T = np.atleast_1d(np.asarray(T, dtype=float))
    P = np.atleast_1d(np.asarray(P, dtype=float))
    T, P = np.broadcast_arrays(T, P)
    result = np.zeros_like(T)
    for i in np.ndindex(T.shape):
        w = IAPWS95(T=float(T[i]), P=float(P[i]))
        result[i] = w.kappa
    return float(result) if result.ndim == 0 else result.squeeze()


if __name__ == "__main__":
    T_test = 298.15
    P_test = 0.1

    print(f"Water properties at T={T_test} K, P={P_test} MPa:")
    print(f"  Backend: {'IAPWS-95' if USE_IAPWS95 else 'IAPWS-IF97 Region 1'}")
    print(f"  rho     = {rho_w(T_test, P_test):.4f} kg/m3")
    print(f"  kappa_T = {kappa_T(T_test, P_test):.6f} 1/MPa")
    print(f"  V1*     = {V1_star(T_test, P_test):.4f} cm3/mol")

    T_res = 373.15
    P_res = 30.0
    print(f"\nWater properties at T={T_res} K, P={P_res} MPa:")
    print(f"  rho     = {rho_w(T_res, P_res):.4f} kg/m3")
    print(f"  kappa_T = {kappa_T(T_res, P_res):.6f} 1/MPa")
    print(f"  V1*     = {V1_star(T_res, P_res):.4f} cm3/mol")
