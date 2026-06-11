/// Named constants for the VLP module, mirroring the constants block at the
/// top of pyrestoolbox/nodal/nodal.py. Keep the two files in sync - the
/// Python module is the authoritative reference.

// ============================================================================
//  Physical constants
// ============================================================================
pub const MW_AIR: f64 = 28.97; // Molecular weight of air (lb/lb-mol)
pub const G_FT: f64 = 32.174; // Gravitational acceleration (ft/s^2)
pub const GC: f64 = 32.174; // Force conversion constant (lbm*ft/(lbf*s^2))
pub const G_SI: f64 = 9.80665; // Gravitational acceleration (m/s^2)
pub const R_GAS: f64 = 10.732; // Gas constant (psia*ft^3/(lb-mol*degR))
pub const P_ATM_PA: f64 = 101325.0; // Atmospheric pressure (Pa)

// ============================================================================
//  Unit conversions
// ============================================================================
pub const PSI_TO_PA: f64 = 6894.757; // psi -> Pa
pub const FT_TO_M: f64 = 0.3048; // ft -> m
pub const IN_TO_M: f64 = 0.0254; // inch -> m
pub const LBFT3_TO_KGM3: f64 = 16.01846; // lb/ft^3 -> kg/m^3
pub const CP_TO_PAS: f64 = 0.001; // cP -> Pa*s
pub const DYNECM_TO_NM: f64 = 0.001; // dyne/cm -> N/m
pub const CP_TO_LBFTS: f64 = 6.7197e-4; // cP -> lb/(ft*s)
pub const LB_TO_KG: f64 = 0.453592; // lb -> kg
pub const DYNCM_PER_LBM_S2: f64 = 453.592; // dyne/cm per lbm/s^2

// ============================================================================
//  Standard fluid and unit constants
// ============================================================================
pub const RHO_AIR_STC: f64 = 0.0765; // Air density at standard conditions (lb/scf)
pub const RHO_FW: f64 = 62.4; // Freshwater density (lb/ft^3)
pub const FT3_PER_BBL: f64 = 5.615; // ft^3 per barrel
pub const SEC_PER_DAY: f64 = 86400.0; // Seconds per day
pub const IN2_PER_FT2: f64 = 144.0; // in^2 per ft^2
pub const FW_GRAD: f64 = 0.433; // Freshwater pressure gradient (psi/ft)

// ============================================================================
//  HB (Hagedorn and Brown 1965) holdup constants
//  Hagedorn, A.R. and Brown, K.E. (1965), JPT 17(4)
// ============================================================================

// Dimensionless group coefficients (also used in BB liquid velocity number)
pub const DG_VEL: f64 = 1.938; // Velocity number coefficient
pub const DG_DIAM: f64 = 120.872; // Pipe diameter number coefficient
pub const DG_VISC: f64 = 0.15726; // Viscosity number coefficient

// CNL polynomial coefficients (5-term)
pub const CNL_C0: f64 = -2.69851;
pub const CNL_C1: f64 = 0.1584095;
pub const CNL_C2: f64 = -0.5509976;
pub const CNL_C3: f64 = 0.5478492;
pub const CNL_C4: f64 = -0.1219458;

// YL/NSI holdup polynomial coefficients (5-term)
pub const YL_C0: f64 = -0.10306578;
pub const YL_C1: f64 = 0.617774;
pub const YL_C2: f64 = -0.632946;
pub const YL_C3: f64 = 0.29598;
pub const YL_C4: f64 = -0.0401;

// S-correction polynomial coefficients (5-term)
pub const SC_C0: f64 = 0.9116257;
pub const SC_C1: f64 = -4.821756;
pub const SC_C2: f64 = 1232.25;
pub const SC_C3: f64 = -22253.58;
pub const SC_C4: f64 = 116174.3;

// S-correction thresholds
pub const SC_F2_CLAMP: f64 = 0.012; // f2 lower clamp
pub const SC_F2_FLOOR: f64 = 0.001; // f2 floor where S = 1.0

// f1 exponents
pub const F1_VG_EXP: f64 = 0.575; // Gas velocity number exponent
pub const F1_P_ATM: f64 = 14.7; // Atmospheric pressure reference (psia)

// f2 exponents
pub const F2_VISC_EXP: f64 = 0.38; // Viscosity number exponent
pub const F2_DIAM_EXP: f64 = 2.14; // Diameter number exponent

// Orkiszewski (1967), JPT 19(6) bubble flow
pub const ORK_VS: f64 = 0.8; // Bubble slip velocity (ft/s)
pub const ORK_LB_A: f64 = 1.071; // Bubble flow boundary coefficient
pub const ORK_LB_B: f64 = -0.2218; // Bubble flow boundary velocity coefficient
pub const ORK_LB_MIN: f64 = 0.13; // Minimum bubble flow boundary

// HB Reynolds/friction dimensional constants
// Hagedorn and Brown (1965), JPT 17(4): two-phase Reynolds number
// NRe = 2.2e-2 * w / (D * muL^HL * muG^(1-HL)) with w in lbm/day.
// Mass flow here is in lbm/s, so the constant is 2.2e-2 * 86400 = 1900.8
pub const HB_RE_K: f64 = 1900.8; // Reynolds number dimensional constant (lbm/s basis)
pub const HB_FRIC_K: f64 = 7.413e10; // Friction gradient dimensional constant

// ============================================================================
//  BB (Beggs and Brill 1973) constants
//  Beggs, H.D. and Brill, J.P. (1973), JPT 25(5)
// ============================================================================

// Flow pattern boundaries (L1-L4): coefficient, exponent
pub const BB_L1_A: f64 = 316.0;
pub const BB_L1_B: f64 = 0.302;
pub const BB_L2_A: f64 = 0.0009252;
pub const BB_L2_B: f64 = -2.4684;
pub const BB_L3_A: f64 = 0.10;
pub const BB_L3_B: f64 = -1.4516;
pub const BB_L4_A: f64 = 0.5;
pub const BB_L4_B: f64 = -6.738;

// Horizontal holdup (a, b, c per regime)
pub const BB_HL_SEG: (f64, f64, f64) = (0.98, 0.4846, 0.0868);
pub const BB_HL_INT: (f64, f64, f64) = (0.845, 0.5351, 0.0173);
pub const BB_HL_DIS: (f64, f64, f64) = (1.065, 0.5824, 0.0609);

// Inclination correction (e, f, g, h per regime)
pub const BB_IC_SEG: (f64, f64, f64, f64) = (0.011, -3.7680, 3.5390, -1.6140);
pub const BB_IC_INT: (f64, f64, f64, f64) = (2.960, 0.3050, -0.4473, 0.0978);

// Payne et al. (1979), JPT 31(9) upward flow correction
pub const BB_PAYNE: f64 = 0.924;

// Friction ratio S-factor polynomial
pub const BB_SF_C0: f64 = -0.0523;
pub const BB_SF_C1: f64 = 3.182;
pub const BB_SF_C2: f64 = -0.8725;
pub const BB_SF_C3: f64 = 0.01853;

// Friction S-factor clamp bounds
pub const BB_SF_LO: f64 = -5.0;
pub const BB_SF_HI: f64 = 5.0;

// ============================================================================
//  Gray (1974) holdup and roughness constants
//  Gray, H.E. (1974), "Vertical Flow Correlation in Gas Wells", User's Manual
//  for API 14B Subsurface Controlled Safety Valve Sizing Computer Program.
//  Gas void fraction: fg = (1 - exp(A)) / (R + 1) with R = vsl/vsg,
//  A = -2.314 * (NV * (1 + 205/ND))^B,
//  B = 0.0814 * (1 - 0.0554 * ln(1 + 730*R/(R+1)))
// ============================================================================

pub const GRAY_A_COEF: f64 = 2.314; // Holdup exponent leading coefficient
pub const GRAY_B1: f64 = 0.0814; // Exponent B coefficient
pub const GRAY_B2: f64 = 0.0554; // Exponent B log coefficient
pub const GRAY_B3: f64 = 730.0; // Exponent B log argument coefficient
pub const GRAY_ND_COEF: f64 = 205.0; // Diameter number coefficient in A

pub const GRAY_ROUGH_K: f64 = 28.5; // Effective roughness coefficient
pub const GRAY_R_THRESH: f64 = 0.007; // R threshold for roughness interpolation
pub const GRAY_ROUGH_FLOOR: f64 = 2.77e-5; // Minimum effective roughness (ft)

// ============================================================================
//  WG (Woldesemayat and Ghajar 2007) drift-flux constants
//  Woldesemayat, M.A. and Ghajar, A.J. (2007), Int. J. Multiphase Flow 33(4)
// ============================================================================

pub const WG_DRIFT_K: f64 = 2.9; // Zuber-Findlay drift velocity coefficient
pub const WG_INCL_K: f64 = 1.22; // Inclination correction coefficient

// WG Chisholm (1967) C values (Lockhart-Martinelli)
pub const WG_CHIS_TT: f64 = 20.0; // Turbulent-turbulent
pub const WG_CHIS_LT: f64 = 12.0; // Laminar-turbulent
pub const WG_CHIS_TL: f64 = 10.0; // Turbulent-laminar
pub const WG_CHIS_LL: f64 = 5.0; // Laminar-laminar
