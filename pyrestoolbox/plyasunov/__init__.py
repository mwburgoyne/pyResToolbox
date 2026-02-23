"""
Plyasunov apparent molar volume model for dissolved gases in water.

Provides V_phi(gas, T_K, P_MPa) and gas_mw(gas) for:
    H2, N2, CH4, CO2, C2H6, C3H8, NC4H10, H2S

Based on:
    Plyasunov, A.V. Fluid Phase Equilibria (2019, 2020, 2021) Parts I-III.
    IAPWS-IF97 Region 1 for pure water properties.
"""

from .plyasunov_model import V_phi, gas_mw, V2_inf, B12, A12_inf, MW_GAS
