"""Named constants for the oil module — extracted from inline magic numbers."""

import numpy as np

# API / SG conversion:  SG = 141.5 / (API + 131.5)
_API_NUMER = 141.5
_API_DENOM = 131.5

# Darcy flow constants
_DARCY_RADIAL = 0.00708   # Darcy radial flow constant (STB·cP / (mD·ft·psi·day))
_DARCY_SKIN_OFFSET = 0.75 # Shape-factor offset in ln(re/rw) + S - 0.75

# Vogel (1968) IPR constants
_VOGEL_AOF_DENOM = 1.8    # qmax = J·Pb / 1.8
_VOGEL_LIN = 0.2          # Linear coefficient
_VOGEL_QUAD = 0.8         # Quadratic coefficient

# Standing (1947) Pb / Rs / Bo — Eqs 1.63-1.72
_STAN_T_COEFF = 0.00091   # Temperature coefficient in exponent
_STAN_API_COEFF = 0.0125  # API coefficient in exponent (sign applied at use)
_STAN_DENOM = 18.2        # Denominator constant
_STAN_SG_EXP = 0.83       # SG exponent on (Rs/sg_g) term
_STAN_OFFSET = 1.4        # Offset constant subtracted inside expression

# Velarde, Blasingame & McCain (1997) — Eq 13-14
_VEL_X_C0 = 0.013098     # x-term coefficient
_VEL_X_T_EXP = 0.282372  # Temperature exponent in x-term
_VEL_X_API_C = 8.2e-6    # API coefficient in x-term
_VEL_X_API_EXP = 2.176124  # API exponent in x-term
_VEL_PB_MULT = 1091.47   # Pb multiplier (psig)
_VEL_RS_EXP = 0.081465   # Rs exponent
_VEL_SG_EXP = -0.161488  # sg_sp exponent
_VEL_PB_OFFSET = 0.740152  # Offset inside Pb brackets
_VEL_PB_POW = 5.354891   # Overall Pb power

# Velarde Rs sub-pressure coefficients (A/B/C) — Eqs 3.8a-3.8f
_VEL_RS_A = [9.73e-7, 1.672608, 0.929870, 0.247235, 1.056052]
_VEL_RS_B = [0.022339, -1.004750, 0.337711, 0.132795, 0.302065]
_VEL_RS_C = [0.725167, -1.485480, -0.164741, -0.091330, 0.047094]

# Valko-McCain (2003) Pb — Eq 2-1
_VALMC_C = [
    [-5.48, 1.27, 4.51, -0.7835],
    [-0.0378, -0.0449, -10.84, 6.23e-3],
    [0.281, 4.36e-4, 8.39, -1.22e-5],
    [-0.0206, -4.76e-6, -2.34, 1.03e-8],
]
_VALMC_LNPB_A0 = 7.475   # Intercept in ln(Pb) polynomial
_VALMC_LNPB_A1 = 0.713   # Linear coefficient
_VALMC_LNPB_A2 = 0.0075  # Quadratic coefficient

# Beggs-Robinson (1975) dead/live oil viscosity — Eq 3.23
_BR_Z0 = 3.0324           # Dead oil Z constant
_BR_Z_API = -0.02023      # Dead oil Z API coefficient
_BR_T_EXP = -1.163        # Dead oil temperature exponent
_BR_A_MULT = 10.715       # Live oil A multiplier
_BR_A_EXP = -0.515        # Live oil A exponent
_BR_A_OFFSET = 100.0      # Live oil A Rs offset
_BR_B_MULT = 5.44         # Live oil B multiplier
_BR_B_EXP = -0.338        # Live oil B exponent
_BR_B_OFFSET = 150.0      # Live oil B Rs offset

# Petrosky-Farshad (1993) undersaturated viscosity — Eq 3.24
_PF_POLY = (-1.0146, 1.3322, -0.4876, -1.15036)  # ln(uob) polynomial coeffs
_PF_P_COEFF = 1.3449e-3   # Pressure coefficient

# Standing-Witte-McCain-Hill (1995) density — Eqs 3.17-3.19
_SWMH_RHOPO_INIT = 52.8   # Initial pseudo-liquid density constant
_SWMH_RHOPO_RS = -0.01    # Rs coefficient in initial estimate (rho = 52.8 - 0.01*Rs)
_SWMH_RHOA = np.array([-49.8930, 85.0149, -3.70373, 0.0479818, 2.98914, -0.0356888])  # Eq 3.18c
_SWMH_MASS_NUMER_OIL = 4600   # Oil mass numerator in Eq 3.18b (scf·SG + 4600·SG_o)
_SWMH_MASS_DENOM = 73.71      # Denominator constant in Eq 3.18b
# sg_g path (no sg_sp) — Eq 3.17e
_SWMH_SG_RHOA_A = 38.52
_SWMH_SG_RHOA_B = -0.00326
_SWMH_SG_RHOA_C = 94.75
_SWMH_SG_RHOA_D = -33.93
# Pressure correction — Eq 3.19d
_SWMH_DP_A = 0.167
_SWMH_DP_B = 16.181
_SWMH_DP_C = -0.0425
_SWMH_DP_D = 0.01    # coefficient on squared term (applied as -0.01)
_SWMH_DP_E = 0.299
_SWMH_DP_F = 263.0
_SWMH_DP_G = -0.0603
# Temperature correction — Eq 3.19f
_SWMH_DT_A = 0.00302
_SWMH_DT_B = 1.505
_SWMH_DT_C = -0.951
_SWMH_DT_D = 0.938
_SWMH_DT_E = 0.0216
_SWMH_DT_F = -0.0233
_SWMH_DT_G = -0.0161
_SWMH_DT_H = 0.475

# cofb McCain (Eq 3.13) — undersaturated oil compressibility
_COFB_C = [
    [3.011, -0.0835, 3.51, 0.327, -1.918, 2.52],
    [-2.6254, -0.259, -0.0289, -0.608, -0.642, -2.73],
    [0.497, 0.382, -0.0584, 0.0911, 0.154, 0.429],
]
_COFB_A0 = 2.434     # Intercept in ln(cofb) polynomial
_COFB_A1 = 0.475     # Linear coefficient
_COFB_A2 = 0.048     # Quadratic coefficient

# sg_evolved_gas — McCain & Hill (1995, SPE 30773)
_SGEVOL_THRESHOLD = 314.7  # Pressure threshold between high/low correlations
_SGEVOL_HIGH = [  # p > 314.7 psia
    0, -208.0797, 22885, -0.000063641, 3.38346,
    -0.000992, -0.000081147, -0.001956, 1.081956, 0.394035,
]
_SGEVOL_LOW = [  # p <= 314.7 psia
    0, -214.0887, 9971, -0.001303, 3.12715,
    -0.001495, -0.000085243, -0.003667, 1.47156, 0.714002,
]

# Bo Standing — Eq 1.69
_BO_STAN_A = 0.972
_BO_STAN_B = 1.47e-4
_BO_STAN_C = 1.25     # degf multiplier inside ()
_BO_STAN_EXP = 1.175  # Overall exponent

# Bo McCain — Eq 3.21
_BO_MC_WDEN = 62.372       # Water density used in McCain Bo (lb/cuft at 60°F)
_BO_MC_RS_COEFF = 0.01357  # coefficient on Rs·sg_g term

# Valko-McCain (2003) stock-tank gas SG — Eq 4-2
_SG_ST_C = [
    [-17.275, -0.3354, 3.705, -155.52, 2.085],
    [7.9597, -0.3346, -0.4273, 629.61, -7.097e-2],
    [-1.1013, 0.1956, 1.818e-2, -957.38, 9.859e-4],
    [2.7735e-2, -3.4374e-2, -3.459e-4, 647.57, -6.312e-6],
    [3.2287e-3, 2.08e-3, 2.505e-6, -163.26, 1.4e-8],
]
_SG_ST_POLY = (1.219, 0.198, 0.0845, 0.03, 0.003)  # sg_st = a + bZ + cZ^2 + dZ^3 + eZ^4

# Valko-McCain (2003) stock-tank Rs — Eq 3-2
_RS_ST_C = [
    [-8.005, 1.224, -1.587],
    [2.7, -0.5, 0.0441],
    [-0.161, 0, -2.29e-5],
]
_RS_ST_POLY = (3.955, 0.83, -0.024, 0.075)  # rs_st = a + bZ + cZ^2 + dZ^3
