"""
Plyasunov A12-infinity model for computing V2-infinity (infinite-dilution
partial molar volume) of dissolved gases in water.

Covers 8 gases using unified equation form from:
  - Part II (Plyasunov 2020): CO2, C2H6, C3H8, n-C4H10
  - Part III (Plyasunov 2021): H2S
  - Part IV (Plyasunov 2021): H2, N2, CH4

Master equation:
    A12_inf = 1 + rho1* * [a0 + a1*rho1* + a2*(rho1*)^2 + ... + a5*(rho1*)^5]

Recovery of V2_inf:
    V2_inf(T,P) = A12_inf(T, rho1*) * kappa_T(T,P) * R * T

Units:
    T in K, P in MPa, rho1* in kg/m3, V2_inf in cm3/mol, kappa_T in 1/MPa

References:
    Plyasunov, A.V. (2019). "Values of the apparent molar volumes V_phi and
    the osmotic coefficients phi of NaCl(aq) at infinite dilution...Part I:
    Aqueous solutions of H2, N2, and CH4." Fluid Phase Equilibria, 496, 43-51.

    Plyasunov, A.V. (2020). "...Part II: CO2, C2H6, C3H8, n-C4H10."
    Fluid Phase Equilibria, 523, 112757.

    Plyasunov, A.V. (2021). "...Part III: H2S."
    Fluid Phase Equilibria, 530, 112883.
"""

import numpy as np
from .water_properties import rho_w, kappa_T, TC_WATER, MW_WATER

# Constants
R_CM3_MPA = 8.314462   # cm3*MPa/(mol*K)
NA = 6.02214076e23     # Avogadro's number
OMEGA = 1e-3 / MW_WATER  # converts B12 from cm3/mol to m3/kg


# ============================================================================
# B12 CROSS SECOND VIRIAL COEFFICIENTS
# ============================================================================

def _B12_polynomial(T, coeffs):
    """
    B12 for simple fluids (H2, N2, CH4) from Part IV Eq. 5:
    B12 = sum( ai * (T/100)^bi )  [cm3/mol]
    coeffs: list of (a, b) tuples, 4 terms
    """
    T100 = T / 100.0
    result = 0.0
    for a, b in coeffs:
        result += a * T100**b
    return result


def _B12_CO2(T):
    """
    B12 for CO2 from Part II Eq. 8 (Hellmann 2019):
    B12 = b1 + b2/(T*)^0.5 + b3/T* + b4/(T*)^3 + b5/(T*)^6 + b6/(T*)^(21/2)
    where T* = T/100, result in cm3/mol
    """
    Tstar = T / 100.0
    b = [15.210, 149.72, -534.54, -2234.6, -13017.0, -39482.0]
    return (b[0]
            + b[1] / Tstar**0.5
            + b[2] / Tstar
            + b[3] / Tstar**3
            + b[4] / Tstar**6
            + b[5] / Tstar**10.5)


def _B12_square_well(T, sigma, eps_kB, lam):
    """
    B12 from square-well potential (Part II Eq. 9):
    B12 = (2/3)*pi*NA*sigma12^3 * {1 - (lambda12^3-1)*[exp(epsilon12/(kB*T)) - 1]}

    Parameters:
        sigma: collision diameter in Angstrom
        eps_kB: well depth / kB in K
        lam: well width in collision diameters

    Returns: B12 in cm3/mol
    """
    sigma_cm = sigma * 1e-8
    b_hs = (2.0 / 3.0) * np.pi * NA * sigma_cm**3
    exp_term = np.exp(eps_kB / T) - 1.0
    return b_hs * (1.0 - (lam**3 - 1.0) * exp_term)


def _B12_group_contribution(T, groups):
    """
    B12 for hydrocarbons via group contributions:
    B12(compound) = sum( n_i * B12(group_i) )
    groups: list of (n_copies, sigma, eps_kB, lambda) tuples
    """
    total = 0.0
    for n, sigma, eps_kB, lam in groups:
        total += n * _B12_square_well(T, sigma, eps_kB, lam)
    return total


# Square-well parameters for functional groups (Part II, Table 2)
_SW_CH3 = (2.788, 283.8, 1.430)
_SW_CH2 = (2.226, 271.4, 1.430)
_SW_H2S = (2.85, 650.0, 1.324)  # Part III

# B12 polynomial coefficients (Part IV, Table 1): list of (a, b) tuples
_B12_POLY_H2 = [(33.047, -0.21), (-250.41, -1.50), (285.42, -2.26), (-186.78, -3.21)]
_B12_POLY_N2 = [(156.679, -0.33), (-183.541, -0.57), (-194.330, -1.47), (-154.815, -3.66)]
_B12_POLY_CH4 = [(109.22, -0.2), (-202.52, -0.6), (-235.86, -2.0), (-297.63, -3.0)]


def B12(gas, T):
    """
    Cross second virial coefficient B12 for gas-water interaction.

    Parameters:
        gas: gas name string (case-insensitive)
        T: temperature in K

    Returns:
        B12 in cm3/mol
    """
    gas = gas.upper()
    if gas == 'H2':
        return _B12_polynomial(T, _B12_POLY_H2)
    elif gas == 'N2':
        return _B12_polynomial(T, _B12_POLY_N2)
    elif gas == 'CH4':
        return _B12_polynomial(T, _B12_POLY_CH4)
    elif gas == 'CO2':
        return _B12_CO2(T)
    elif gas == 'C2H6':
        return _B12_group_contribution(T, [(2, *_SW_CH3)])
    elif gas == 'C3H8':
        return _B12_group_contribution(T, [(2, *_SW_CH3), (1, *_SW_CH2)])
    elif gas in ('NC4H10', 'N-C4H10', 'NC4'):
        return _B12_group_contribution(T, [(2, *_SW_CH3), (2, *_SW_CH2)])
    elif gas == 'H2S':
        return _B12_square_well(T, *_SW_H2S)
    else:
        raise ValueError(f"Unknown gas: {gas}")


# ============================================================================
# p_in COEFFICIENT TABLES
# ============================================================================

# Format: p_in[i][n] where i=1..5, n=0..6
# Each gas has a 5x7 array of coefficients

# Paper 8 (Part IV, Table 2) - full precision - supersedes Paper 5
P_IN_H2 = [
    [2.1202801e-3, -6.037483e-3, 1.0195766e-2, -1.0693934e-2, 5.4991895e-3, -1.0878264e-3, 1.1746219e-6],
    [-1.1046217e-5, 3.1667503e-5, -5.3803531e-5, 5.6835389e-5, -2.9567299e-5, 5.9819654e-6, -2.019032e-8],
    [1.9307298e-8, -5.4979556e-8, 9.2491093e-8, -9.7577108e-8, 5.1231461e-8, -1.0724676e-8, 1.2613281e-10],
    [-1.3252458e-11, 3.7335356e-11, -6.1484868e-11, 6.403194e-11, -3.3547995e-11, 7.4031616e-12, -3.4170589e-13],
    [3.011170e-15, -8.3139044e-15, 1.3052306e-14, -1.2814023e-14, 6.0470404e-15, -9.8662924e-16, -4.1237636e-17],
]

P_IN_N2 = [
    [-2.2070753e-3, 6.6067485e-3, -1.2027859e-2, 1.3728692e-2, -7.7929165e-3, 1.7076684e-3, -2.8424026e-6],
    [1.1431298e-5, -3.3833752e-5, 6.1522458e-5, -7.1045896e-5, 4.1165829e-5, -9.3369807e-6, 4.7959089e-8],
    [-1.8029853e-8, 5.3340827e-8, -9.7754731e-8, 1.1514188e-7, -6.8971366e-8, 1.667536e-8, -3.0103987e-10],
    [1.0436147e-11, -3.0813479e-11, 5.7057178e-11, -6.8981444e-11, 4.3377118e-11, -1.2090314e-11, 9.5828102e-13],
    [-1.6676726e-15, 4.8355967e-15, -8.7551827e-15, 1.0376018e-14, -6.124985e-15, 1.3543054e-15, -1.6688106e-18],
]

P_IN_CH4 = [
    [-1.3023759e-3, 3.8481069e-3, -6.8502936e-3, 7.7002273e-3, -4.3454809e-3, 9.5403855e-4, -1.8741835e-6],
    [7.3809403e-6, -2.1303297e-5, 3.7100486e-5, -4.1343101e-5, 2.3418494e-5, -5.2627676e-6, 3.0296735e-8],
    [-1.2366446e-8, 3.5518855e-8, -6.1394616e-8, 6.8386294e-8, -3.9321302e-8, 9.3024856e-9, -1.7821599e-10],
    [7.5953065e-12, -2.1583027e-11, 3.6743322e-11, -4.050275e-11, 2.3523613e-11, -6.1571572e-12, 4.6261515e-13],
    [-1.3716224e-15, 3.7669953e-15, -5.924387e-15, 5.8112619e-15, -2.7435547e-15, 3.8223834e-16, 5.3380661e-17],
]

# Paper 6 (Part II, Table 5) - full precision
P_IN_CO2 = [
    [9.1583134e-4, -2.7160016e-3, 5.0839423e-3, -5.8354599e-3, 3.2543133e-3, -7.0687955e-4, 1.4494891e-6],
    [-3.6106209e-6, 1.1315777e-5, -2.2758954e-5, 2.7797001e-5, -1.6473215e-5, 3.8439408e-6, -2.557928e-8],
    [5.0157576e-9, -1.5637417e-8, 3.200829e-8, -4.1225984e-8, 2.6278766e-8, -6.8414637e-9, 1.6610365e-10],
    [-2.0446824e-12, 6.2444294e-12, -1.3130671e-11, 1.9003321e-11, -1.3986032e-11, 4.6916137e-12, -5.1705069e-13],
    [-1.1881416e-16, 4.4574039e-16, -7.9103055e-16, 2.5194274e-16, 2.4022155e-16, -9.7895994e-17, -1.4094429e-17],
]

P_IN_C2H6 = [
    [-5.1475443e-4, 1.729277e-3, -3.6536764e-3, 4.6194269e-3, -2.7977608e-3, 6.2657963e-4, -7.9712465e-8],
    [3.5428879e-6, -1.0888811e-5, 2.1555401e-5, -2.6562891e-5, 1.5872915e-5, -3.5330389e-6, 1.9105807e-9],
    [-7.076262e-9, 2.150819e-8, -4.1847563e-8, 5.0828275e-8, -3.0091634e-8, 6.7020044e-9, -1.7913067e-11],
    [5.7671272e-12, -1.7389573e-11, 3.3600113e-11, -4.0403444e-11, 2.3837864e-11, -5.457811e-12, 1.0241543e-13],
    [-1.6279128e-15, 4.8764202e-15, -9.3154249e-15, 1.1023937e-14, -6.3736955e-15, 1.4100215e-15, -1.4735524e-17],
]

P_IN_C3H8 = [
    [2.2892124e-3, -6.4740679e-3, 1.1030983e-2, -1.1944689e-2, 6.4760838e-3, -1.3894223e-3, 3.3011188e-6],
    [-1.029341e-5, 3.0090181e-5, -5.3316193e-5, 5.9603698e-5, -3.35325e-5, 7.5922124e-6, -5.4814576e-8],
    [1.5037316e-8, -4.4050725e-8, 7.8975079e-8, -9.0810405e-8, 5.3609476e-8, -1.3301616e-8, 3.3278243e-10],
    [-7.3401228e-12, 2.1306209e-11, -3.7787336e-11, 4.4705147e-11, -2.8278765e-11, 8.5764552e-12, -9.1713991e-13],
    [5.5117215e-16, -1.4000747e-15, 1.7856523e-15, -1.6585346e-15, 6.3567082e-16, 7.7519099e-17, -7.8311547e-17],
]

P_IN_NC4H10 = [
    [2.5603794e-3, -6.9646157e-3, 1.0987245e-2, -1.1110696e-2, 5.7384581e-3, -1.2233119e-3, 4.9218547e-6],
    [-8.9460753e-6, 2.5173896e-5, -4.1659657e-5, 4.4713209e-5, -2.5204428e-5, 6.1085451e-6, -8.1017199e-8],
    [7.0204532e-9, -1.9465263e-8, 3.2883394e-8, -4.0132208e-8, 2.8071928e-8, -9.1449293e-9, 4.8274281e-10],
    [3.6497249e-12, -1.1330901e-11, 2.0340955e-11, -1.7139410e-11, 2.6365251e-12, 3.4703436e-12, -1.2475357e-12],
    [-3.7951579e-15, 1.1392044e-14, -2.075511e-14, 2.2507396e-14, -1.20847e-14, 2.7449667e-15, -1.3968834e-16],
]

# Paper 7 (Part III) - full precision from Table 5
P_IN_H2S = [
    [6.7014672e-4, -1.9403321e-3, 3.5360019e-3, -4.0122786e-3, 2.2524974e-3, -4.8843038e-4, 1.6044152e-6],
    [-3.1805655e-6, 9.5800671e-6, -1.7994973e-5, 2.1087746e-5, -1.2352096e-5, 2.8686985e-6, -2.7620761e-8],
    [5.0089935e-9, -1.4885595e-8, 2.7733075e-8, -3.3416171e-8, 2.0773681e-8, -5.4162153e-9, 1.7656758e-10],
    [-2.5739894e-12, 7.422722e-12, -1.3339226e-11, 1.677997e-11, -1.1669169e-11, 4.0370617e-12, -5.6862157e-13],
    [1.6189502e-16, -3.2073034e-16, 1.4144345e-16, -7.3004629e-17, 5.7935525e-17, -3.11107e-18, 1.3146726e-18],
]

# Map gas names to p_in tables
P_IN_MAP = {
    'H2': P_IN_H2,
    'N2': P_IN_N2,
    'CH4': P_IN_CH4,
    'CO2': P_IN_CO2,
    'C2H6': P_IN_C2H6,
    'C3H8': P_IN_C3H8,
    'NC4H10': P_IN_NC4H10,
    'N-C4H10': P_IN_NC4H10,
    'NC4': P_IN_NC4H10,
    'H2S': P_IN_H2S,
}

# Molecular weights of dissolved gases (g/mol)
MW_GAS = {
    'H2': 2.01588,
    'N2': 28.0134,
    'CH4': 16.0425,
    'CO2': 44.0095,
    'C2H6': 30.069,
    'C3H8': 44.096,
    'NC4H10': 58.122,
    'N-C4H10': 58.122,
    'NC4': 58.122,
    'H2S': 34.081,
}


# ============================================================================
# CORE MODEL
# ============================================================================

def _compute_ai(p_in_row, T, i_index):
    """
    Compute temperature-dependent coefficient a_i from p_in coefficients.

    Eq. 34: ai = pi0*theta^(-0.5) + pi1*theta^(-1) + pi2*theta^(-2) + ...
                 + pi5*theta^(-5) + pi6*theta^(-n)
    where theta = T/Tc, n = 9(i=1), 8(i=2), 7(i=3), 6(i=4,5)
    i_index: 1-based index (1-5)
    """
    theta = T / TC_WATER
    exponents = [-0.5, -1.0, -2.0, -3.0, -4.0, -5.0]
    n_map = {1: 9, 2: 8, 3: 7, 4: 6, 5: 6}
    last_exp = -n_map[i_index]

    result = 0.0
    for k in range(6):
        result += p_in_row[k] * theta**exponents[k]
    result += p_in_row[6] * theta**last_exp

    return result


def A12_inf(gas, T, rho1_star):
    """
    Compute A12_inf = V2_inf/(kappa_T*R*T) at given T and water density.

    Eq. 32: A12_inf = 1 + rho1*[a0 + a1*rho1* + ... + a5*(rho1*)^5]
    where a0 = 2*Omega*B12

    Parameters:
        gas: gas name string
        T: temperature in K
        rho1_star: pure water density in kg/m3

    Returns:
        A12_inf (dimensionless)
    """
    gas = gas.upper()
    if gas not in P_IN_MAP:
        raise ValueError(f"Unknown gas: {gas}. Available: {list(P_IN_MAP.keys())}")

    p_in = P_IN_MAP[gas]

    a0 = 2.0 * OMEGA * B12(gas, T)
    a_coeffs = [a0]
    for i in range(5):
        a_coeffs.append(_compute_ai(p_in[i], T, i + 1))

    rho = rho1_star
    poly_sum = 0.0
    rho_power = 1.0
    for k in range(6):
        poly_sum += a_coeffs[k] * rho_power
        rho_power *= rho

    return 1.0 + rho * poly_sum


def V2_inf(gas, T, P):
    """
    Compute V2_inf (infinite-dilution partial molar volume) at given T, P.

    V2_inf = A12_inf * kappa_T * R * T

    Parameters:
        gas: gas name string (case-insensitive)
        T: temperature in K
        P: pressure in MPa

    Returns:
        V2_inf in cm3/mol
    """
    rho = rho_w(T, P)
    kt = kappa_T(T, P)
    a12 = A12_inf(gas, T, rho)
    return a12 * kt * R_CM3_MPA * T


def V_phi(gas, T, P):
    """
    Apparent molar volume V_phi at infinite dilution.
    Alias for V2_inf.

    Parameters:
        gas: gas name string (case-insensitive).
             Supported: H2, N2, CH4, CO2, C2H6, C3H8, NC4H10, H2S
        T: temperature in K
        P: pressure in MPa

    Returns:
        V_phi in cm3/mol
    """
    return V2_inf(gas, T, P)


def gas_mw(gas):
    """Return molecular weight of gas in g/mol."""
    gas = gas.upper()
    if gas not in MW_GAS:
        raise ValueError(f"Unknown gas: {gas}")
    return MW_GAS[gas]
