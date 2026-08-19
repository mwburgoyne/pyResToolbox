/// Binary Interaction Parameter (BIP) calculations for the S&W VLE engine.
///
/// Provides kij_AQ (aqueous phase) and kij_NA (non-aqueous phase) BIP values
/// for gas-water pairs, plus gas-gas BIPs from literature.

use crate::vle::components::*;

// =============================================================================
// kij_AQ Correlations (Aqueous Phase) — proposed framework (MC-3 alpha,
// freshwater-only kij). Salinity via Sechenov/gamma-phi.
// =============================================================================

/// S&W 1992 Equations 11-12 with errata — general hydrocarbon form.
fn kij_aq_hydrocarbon(t_k: f64, omega: f64, tc: f64, salinity_molal: f64) -> f64 {
    let tr = t_k / tc;
    let cs = salinity_molal;

    // Corrected constants from S&W errata
    let a0 = 1.1120 - 1.7369 * omega.powf(-0.1);
    let a1 = 1.1001 + 0.8360 * omega;
    let a2 = -0.15742 - 1.0988 * omega;

    // Corrected salinity coefficients from errata
    let alpha0 = 0.017407;
    let alpha1 = 0.033516;
    let alpha2 = 0.011478;

    a0 * (1.0 + alpha0 * cs)
        + a1 * tr * (1.0 + alpha1 * cs)
        + a2 * tr * tr * (1.0 + alpha2 * cs)
}

/// CH4-water kij_AQ (this work, rational form, MC-3 alpha).
/// kij = (A + Tr) / (B + C*Tr) where Tr = T/190.60
fn kij_aq_ch4(t_k: f64) -> f64 {
    let tr = t_k / 190.60;
    let (a, b, c) = (-2.1642, 1.7325, 0.2105);
    (a + tr) / (b + c * tr)
}

/// CO2-water kij_AQ (proposed, cubic in T(K), MC-3 alpha).
fn kij_aq_co2_proposed(t_k: f64) -> f64 {
    -1.5127 + 9.4980e-3 * t_k - 2.1680e-5 * t_k * t_k + 1.8500e-8 * t_k * t_k * t_k
}

/// H2S-water kij_AQ (proposed, exp form, MC-3 alpha).
fn kij_aq_h2s_proposed(t_k: f64) -> f64 {
    -70.8170 / t_k + 1540.9516 * (-4532.56 / t_k).exp() + 0.21517
}

/// N2-water kij_AQ (proposed, linear in T(K), MC-3 alpha).
fn kij_aq_n2_proposed(t_k: f64) -> f64 {
    -1.6669 + 3.447873e-3 * t_k
}

/// H2-water kij_AQ (proposed, rational form, MC-3 alpha).
fn kij_aq_h2_proposed(t_k: f64) -> f64 {
    let tr = t_k / 33.145;
    (-14.6157 + tr) / (3.5494 + 0.2230 * tr)
}

/// C2H6-water kij_AQ (proposed, rational, MC-3 alpha).
fn kij_aq_c2h6_proposed(t_k: f64) -> f64 {
    let tr = t_k / 305.40;
    (-1.2683 + tr) / (0.3352 + 1.3035 * tr)
}

/// C3H8-water kij_AQ (proposed, rational, MC-3 alpha).
fn kij_aq_c3h8_proposed(t_k: f64) -> f64 {
    let tr = t_k / 369.80;
    (-1.1491 + tr) / (0.6126 + 1.3183 * tr)
}

/// Freshwater wrapper for HCs without fitted proposed correlations.
/// Uses S&W Eqs 11-12 with cs=0.
fn kij_aq_hc_proposed(comp_idx: usize, t_k: f64) -> f64 {
    let c = &COMPONENT_DB[comp_idx];
    kij_aq_hydrocarbon(t_k, c.omega, c.tc, 0.0)
}

/// Get kij_AQ for any supported gas (proposed framework, freshwater only).
/// Salinity is handled via gamma-phi (Sechenov) externally.
pub fn get_kij_aq_proposed(comp_idx: usize, t_k: f64) -> f64 {
    match comp_idx {
        IDX_H2 => kij_aq_h2_proposed(t_k),
        IDX_CO2 => kij_aq_co2_proposed(t_k),
        IDX_N2 => kij_aq_n2_proposed(t_k),
        IDX_H2S => kij_aq_h2s_proposed(t_k),
        IDX_CH4 => kij_aq_ch4(t_k),
        IDX_C2H6 => kij_aq_c2h6_proposed(t_k),
        IDX_C3H8 => kij_aq_c3h8_proposed(t_k),
        IDX_IC4H10 | IDX_NC4H10 | IDX_IC5H12 | IDX_NC5H12
        | IDX_NC6H14 | IDX_NC7H16 | IDX_NC8H18 | IDX_NC10H22 => {
            kij_aq_hc_proposed(comp_idx, t_k)
        }
        _ => 0.0, // H2O-H2O or unknown
    }
}

// =============================================================================
// Default framework kij_AQ (S&W water alpha, refitted freshwater correlations)
// =============================================================================
// Burgoyne & Nielsen (2026), the published and recommended framework. Salinity
// is NOT in these correlations: it enters through the additive embedded
// delta_kij below. Mirrors KIJ_AQ_DEFAULT in _lib_vle_engine.py.

fn kij_aq_co2_default(t_k: f64) -> f64 {
    -1.5893 + 9.8873e-3 * t_k - 2.2188e-5 * t_k * t_k + 1.8499e-8 * t_k * t_k * t_k
}

fn kij_aq_h2s_default(t_k: f64) -> f64 {
    -74.6914 / t_k + 1348.9615 * (-4504.96 / t_k).exp() + 0.22598
}

fn kij_aq_n2_default(t_k: f64) -> f64 {
    -1.6689 + 3.441589e-3 * t_k
}

fn kij_aq_h2_default(t_k: f64) -> f64 {
    let tr = t_k / 33.145;
    (-14.9412 + tr) / (2.2832 + 0.3893 * tr)
}

fn kij_aq_ch4_default(t_k: f64) -> f64 {
    let tr = t_k / 190.60;
    (-2.1756 + tr) / (1.0388 + 0.6436 * tr)
}

fn kij_aq_c2h6_default(t_k: f64) -> f64 {
    let tr = t_k / 305.40;
    (-1.2669 + tr) / (0.1526 + 1.4335 * tr)
}

fn kij_aq_c3h8_default(t_k: f64) -> f64 {
    let tr = t_k / 369.80;
    (-1.1460 + tr) / (0.5760 + 1.3107 * tr)
}

/// Get kij_AQ for any supported gas (default framework, freshwater only).
pub fn get_kij_aq_default(comp_idx: usize, t_k: f64) -> f64 {
    match comp_idx {
        IDX_H2 => kij_aq_h2_default(t_k),
        IDX_CO2 => kij_aq_co2_default(t_k),
        IDX_N2 => kij_aq_n2_default(t_k),
        IDX_H2S => kij_aq_h2s_default(t_k),
        IDX_CH4 => kij_aq_ch4_default(t_k),
        IDX_C2H6 => kij_aq_c2h6_default(t_k),
        IDX_C3H8 => kij_aq_c3h8_default(t_k),
        IDX_IC4H10 | IDX_NC4H10 | IDX_IC5H12 | IDX_NC5H12
        | IDX_NC6H14 | IDX_NC7H16 | IDX_NC8H18 | IDX_NC10H22 => {
            kij_aq_hc_proposed(comp_idx, t_k)
        }
        _ => 0.0, // H2O-H2O or unknown
    }
}

// =============================================================================
// Embedded salinity delta_kij (paper Table 7)
// =============================================================================
// delta_kij = (a0 + a1*Tr + a2*Tr^2)*m, plus (b0 + b1*Tr)*m^2 for CO2 only.
// Mirrors EMBEDDED_SALINITY_PARAMS_DEFAULT and calc_embedded_delta_kij in
// _lib_vle_engine.py.

/// (tc, a0, a1, a2, b0, b1). b0/b1 are zero for the linear-in-m gases.
fn embedded_delta_params(comp_idx: usize) -> Option<(f64, f64, f64, f64, f64, f64)> {
    match comp_idx {
        IDX_CO2 => Some((304.20, 0.0409, -0.0807, 0.0526, 0.0079, -0.0085)),
        IDX_H2S => Some((373.20, 0.0341, -0.0655, 0.0376, 0.0, 0.0)),
        IDX_CH4 => Some((190.60, 0.1304, -0.1295, 0.0394, 0.0, 0.0)),
        IDX_N2 => Some((126.10, 0.2173, -0.1468, 0.0302, 0.0, 0.0)),
        IDX_H2 => Some((33.145, 0.3658, -0.0625, 0.0030, 0.0, 0.0)),
        IDX_C2H6 => Some((305.40, 0.0813, -0.1287, 0.0646, 0.0, 0.0)),
        IDX_C3H8 => Some((369.80, 0.0606, -0.1165, 0.0772, 0.0, 0.0)),
        IDX_NC4H10 => Some((425.20, 0.0488, -0.1072, 0.0836, 0.0, 0.0)),
        _ => None,
    }
}

/// Embedded salinity correction on kij_AQ. Zero for m <= 0 or an unfitted gas.
pub fn calc_embedded_delta_kij(comp_idx: usize, t_k: f64, salinity_molal: f64) -> f64 {
    if salinity_molal <= 0.0 {
        return 0.0;
    }
    match embedded_delta_params(comp_idx) {
        None => 0.0,
        Some((tc, a0, a1, a2, b0, b1)) => {
            let tr = t_k / tc;
            let m = salinity_molal;
            (a0 + a1 * tr + a2 * tr * tr) * m + (b0 + b1 * tr) * m * m
        }
    }
}

/// Which BIP set and water alpha the flash should use.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Framework {
    /// S&W water alpha + refitted freshwater kij + embedded delta_kij.
    Default,
    /// Mathias-Copeman 3-parameter water alpha + freshwater kij + gamma-phi.
    Mc3,
}

// =============================================================================
// kij_NA Values (Non-Aqueous Phase) — S&W 1992 Table 5 + this work
// =============================================================================

/// Constant kij_NA values. Index matches COMPONENT_NAMES.
/// H2O entry is 0.0 (unused).
const KIJ_NA_TABLE: [f64; NUM_COMPONENTS] = [
    0.0,    // H2O (placeholder)
    0.468,  // H2 (this work)
    0.1896, // CO2
    0.1610, // H2S (this work, constant; replaces S&W Eq 17)
    0.4778, // N2
    0.4850, // CH4
    0.4920, // C2H6
    0.5525, // C3H8 (S&W 1992 Table 5; 0.5070 was a transcription error)
    0.5091, // iC4H10 (not in S&W Table 5; set equal to nC4H10)
    0.5091, // nC4H10 (S&W 1992 Table 5; 0.5080 was a transcription error)
    0.5000, // iC5H12 (Whitson & Brule 2000 Table 9.3, C5+)
    0.5000, // nC5H12 (Whitson & Brule 2000 Table 9.3, C5+)
    0.5000, // nC6H14 (Whitson & Brule 2000 Table 9.3, C5+)
    0.5000, // nC7H16 (Whitson & Brule 2000 Table 9.3, C5+)
    0.5000, // nC8H18 (Whitson & Brule 2000 Table 9.3, C5+)
    0.5000, // nC10H22 (Whitson & Brule 2000 Table 9.3, C5+)
];

/// Get kij_NA for a given component index.
pub fn get_kij_na(comp_idx: usize) -> f64 {
    if comp_idx < NUM_COMPONENTS {
        KIJ_NA_TABLE[comp_idx]
    } else {
        0.0
    }
}

// =============================================================================
// Gas-Gas BIPs (Literature Values)
// =============================================================================

/// Gas-gas BIP lookup. Returns 0.0 for unknown or self-pairs.
pub fn get_gas_gas_bip(idx_a: usize, idx_b: usize) -> f64 {
    if idx_a == idx_b {
        return 0.0;
    }
    // Normalize order: smaller index first
    let (a, b) = if idx_a < idx_b {
        (idx_a, idx_b)
    } else {
        (idx_b, idx_a)
    };

    // Lookup table indexed by (smaller_idx, larger_idx)
    // Only gas-gas pairs (not water)
    match (a, b) {
        (IDX_H2, IDX_CO2) => 0.0,
        (IDX_H2, IDX_H2S) => 0.0,
        (IDX_H2, IDX_N2) => 0.0,
        (IDX_H2, IDX_CH4) => 0.0,
        (IDX_H2, IDX_C2H6) => 0.0,
        (IDX_H2, IDX_C3H8) => 0.0,
        (IDX_H2, IDX_IC4H10) => 0.0,
        (IDX_H2, IDX_NC4H10) => 0.0,
        (IDX_CO2, IDX_H2S) => 0.097,
        (IDX_CO2, IDX_N2) => -0.02,
        (IDX_CO2, IDX_CH4) => 0.12,
        (IDX_CO2, IDX_C2H6) => 0.13,
        (IDX_CO2, IDX_C3H8) => 0.135,
        (IDX_CO2, IDX_IC4H10) => 0.13,
        (IDX_CO2, IDX_NC4H10) => 0.13,
        (IDX_H2S, IDX_N2) => 0.17,
        (IDX_H2S, IDX_CH4) => 0.08,
        (IDX_H2S, IDX_C2H6) => 0.085,
        (IDX_H2S, IDX_C3H8) => 0.08,
        (IDX_N2, IDX_CH4) => 0.036,
        (IDX_N2, IDX_C2H6) => 0.04,
        (IDX_N2, IDX_C3H8) => 0.08,
        (IDX_CH4, IDX_C2H6) => 0.0026,
        (IDX_CH4, IDX_C3H8) => 0.014,
        (IDX_CH4, IDX_IC4H10) => 0.02,
        (IDX_CH4, IDX_NC4H10) => 0.02,
        (IDX_C2H6, IDX_C3H8) => 0.001,
        (IDX_C2H6, IDX_NC4H10) => 0.01,
        (IDX_C3H8, IDX_NC4H10) => 0.003,
        _ => 0.0,
    }
}

// =============================================================================
// Sechenov Coefficient (S&W 1992 Equation 8)
// =============================================================================

/// Soreide-Whitson (1992) Equation 8 — Sechenov coefficient correlation.
///
/// ks = 0.13163 + 4.45e-4*Tb - 7.692e-4*T_F + 2.6614e-6*T_F^2 - 2.612e-9*T_F^3
///
/// where T_F is temperature in Fahrenheit and Tb is normal boiling point in K.
///
/// Returns ks on log10 basis (per S&W Eq. 2).
pub fn sw_equation_8_ks(t_c: f64, tb_k: f64) -> f64 {
    let t_f = t_c * 9.0 / 5.0 + 32.0;
    0.13163 + 4.45e-4 * tb_k - 7.692e-4 * t_f
        + 2.6614e-6 * t_f * t_f
        - 2.612e-9 * t_f * t_f * t_f
}

/// Get Sechenov coefficient ks for a given gas (proposed framework).
///
/// For CO2 and H2S in the full proposed framework, specialized models are used
/// on the Python side. This Rust implementation provides the S&W Eq 8 fallback
/// plus the N2 +0.02 offset.
pub fn get_sechenov_ks_sw(comp_idx: usize, t_k: f64) -> f64 {
    if comp_idx >= NUM_COMPONENTS || comp_idx == IDX_H2O {
        return 0.0;
    }
    let t_c = t_k - 273.15;
    let tb_k = COMPONENT_DB[comp_idx].tb;
    let mut ks = sw_equation_8_ks(t_c, tb_k);
    // N2: empirical +0.02 offset (proposed mode)
    if comp_idx == IDX_N2 {
        ks += 0.02;
    }
    ks
}

// =============================================================================
// Full kij matrix construction
// =============================================================================

/// Build N x N kij matrix for a given set of component indices and mode.
///
/// `mode_aq`: true for AQ mode (gas-water BIPs from kij_AQ), false for NA mode.
/// Returns flattened row-major N x N matrix.
pub fn build_kij_matrix(
    comp_indices: &[usize],
    t_k: f64,
    mode_aq: bool,
    framework: Framework,
    salinity_molal: f64,
) -> Vec<f64> {
    let nc = comp_indices.len();
    let mut kij = vec![0.0; nc * nc];

    for i in 0..nc {
        for j in (i + 1)..nc {
            let idx_i = comp_indices[i];
            let idx_j = comp_indices[j];

            let val = if idx_i == IDX_H2O || idx_j == IDX_H2O {
                // Gas-water pair
                let gas_idx = if idx_i == IDX_H2O { idx_j } else { idx_i };
                if mode_aq {
                    match framework {
                        Framework::Mc3 => get_kij_aq_proposed(gas_idx, t_k),
                        // Salinity rides on the BIP in the default framework,
                        // so it is added here rather than through gamma.
                        Framework::Default => {
                            get_kij_aq_default(gas_idx, t_k)
                                + calc_embedded_delta_kij(gas_idx, t_k, salinity_molal)
                        }
                    }
                } else {
                    get_kij_na(gas_idx)
                }
            } else {
                // Gas-gas pair
                get_gas_gas_bip(idx_i, idx_j)
            };

            kij[i * nc + j] = val;
            kij[j * nc + i] = val;
        }
    }
    kij
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kij_aq_ch4() {
        // Verify CH4 kij_AQ at 373.15 K
        let kij = kij_aq_ch4(373.15);
        // Tr = 373.15/190.60 = 1.9580
        // kij = (-2.1642 + 1.9580) / (1.7325 + 0.2105*1.9580)
        assert!((kij - (-0.09566)).abs() < 0.001);
    }

    #[test]
    fn test_kij_na_ch4() {
        assert!((get_kij_na(IDX_CH4) - 0.485).abs() < 1e-10);
    }

    #[test]
    fn test_gas_gas_bip_symmetric() {
        assert_eq!(
            get_gas_gas_bip(IDX_CH4, IDX_CO2),
            get_gas_gas_bip(IDX_CO2, IDX_CH4)
        );
    }

    #[test]
    fn test_sw_eq8_ks() {
        // Test at 25C for H2 (Tb=20.3 K)
        let ks = sw_equation_8_ks(25.0, 20.3);
        assert!(ks > 0.0);
        assert!(ks < 0.5);
    }
}
