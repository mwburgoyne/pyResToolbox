//! Binary Interaction Parameter (BIP) calculations for the S&W VLE engine.
//!
//! Provides kij_AQ (aqueous phase) and kij_NA (non-aqueous phase) BIP values
//! for gas-water pairs, plus gas-gas BIPs from literature.

use crate::vle::components::*;

// =============================================================================
// kij_AQ Correlations (Aqueous Phase), 'default' framework
// =============================================================================
// S&W water alpha with freshwater kij_AQ refitted on procedure 1 (S&W 1992
// Eqs 10a/10b: one equilibrium, kij_AQ in the aqueous phase, kij_NA in the
// vapour carrying its water). Burgoyne & Nielsen (2026) forms and data.
// Salinity is NOT in these correlations: it enters through the additive
// embedded delta_kij below. Mirrors KIJ_AQ_DEFAULT and the kij_aq_*_dropin
// functions in _lib_vle_engine.py. The parameters are strongly correlated:
// keep every digit.

/// S&W 1992 Equations 11-12 with errata, general hydrocarbon form.
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

/// CO2: cubic in T (K). Port of kij_aq_co2_dropin.
fn kij_aq_co2_default(t_k: f64) -> f64 {
    -1.5502105835458762 + 0.009322719516912797 * t_k
        + -1.977702057788579e-05 * t_k.powf(2.0)
        + 1.527732067553674e-08 * t_k.powf(3.0)
}

/// H2S: a/T + b*exp(c/T) + d. Port of kij_aq_h2s_dropin.
fn kij_aq_h2s_default(t_k: f64) -> f64 {
    -82.78098013402487 / t_k + 1739490.5143550164 * (-8585.063329478664 / t_k).exp()
        + 0.2518306228082515
}

/// N2: linear in T (K). Port of kij_aq_n2_dropin.
fn kij_aq_n2_default(t_k: f64) -> f64 {
    -1.6210963967981318 + 0.003295272503921554 * t_k
}

/// H2: rational in Tr = T/33.145. Port of kij_aq_h2_dropin.
fn kij_aq_h2_default(t_k: f64) -> f64 {
    let tr = t_k / 33.145;
    (-16.019825767516267 + tr) / (0.7312362834651882 + 0.6813619680172914 * tr)
}

/// He: rational in Tr = T/5.1953 (NIST Tc). Port of kij_aq_he_dropin.
fn kij_aq_he_default(t_k: f64) -> f64 {
    let tr = t_k / 5.1953;
    (-71.61550493445799 + tr) / (-0.1644070191372901 + 0.5017371389503704 * tr)
}

/// CH4: rational in Tr = T/190.60. Port of kij_aq_ch4_dropin.
fn kij_aq_ch4_default(t_k: f64) -> f64 {
    let tr = t_k / 190.60;
    (-2.1976782374118784 + tr) / (0.7150153669761119 + 0.8945029639230546 * tr)
}

/// C2H6: rational in Tr = T/305.40. Port of kij_aq_c2h6_dropin.
fn kij_aq_c2h6_default(t_k: f64) -> f64 {
    let tr = t_k / 305.40;
    (-1.271347888293827 + tr) / (-0.02927136624699385 + 1.6371834535180028 * tr)
}

/// C3H8: rational in Tr = T/369.80. Port of kij_aq_c3h8_dropin.
fn kij_aq_c3h8_default(t_k: f64) -> f64 {
    let tr = t_k / 369.80;
    (-1.1567778736282177 + tr) / (0.2089459674695665 + 1.81066287107836 * tr)
}

/// Freshwater S&W Eqs 11-12 (cs = 0) for butanes and heavier.
fn kij_aq_hc_freshwater(comp_idx: usize, t_k: f64) -> f64 {
    let c = &COMPONENT_DB[comp_idx];
    kij_aq_hydrocarbon(t_k, c.omega, c.tc, 0.0)
}

/// Freshwater kij_AQ for any supported gas ('default' framework).
pub fn get_kij_aq_default(comp_idx: usize, t_k: f64) -> f64 {
    match comp_idx {
        IDX_H2 => kij_aq_h2_default(t_k),
        IDX_CO2 => kij_aq_co2_default(t_k),
        IDX_N2 => kij_aq_n2_default(t_k),
        IDX_H2S => kij_aq_h2s_default(t_k),
        IDX_CH4 => kij_aq_ch4_default(t_k),
        IDX_C2H6 => kij_aq_c2h6_default(t_k),
        IDX_C3H8 => kij_aq_c3h8_default(t_k),
        IDX_HE => kij_aq_he_default(t_k),
        IDX_IC4H10 | IDX_NC4H10 | IDX_IC5H12 | IDX_NC5H12
        | IDX_NC6H14 | IDX_NC7H16 | IDX_NC8H18 | IDX_NC10H22 => {
            kij_aq_hc_freshwater(comp_idx, t_k)
        }
        _ => 0.0, // H2O-H2O or unknown
    }
}

// =============================================================================
// Embedded salinity delta_kij (procedure-1 refit)
// =============================================================================
// delta_kij = (a0 + a1*Tr + a2*Tr^2)*m, plus (b0 + b1*Tr)*m^2 for CO2 only.
// Mirrors EMBEDDED_SALINITY_PARAMS_DEFAULT and calc_embedded_delta_kij in
// _lib_vle_engine.py.

/// (tc, a0, a1, a2, b0, b1). b0/b1 are zero for the linear-in-m gases.
fn embedded_delta_params(comp_idx: usize) -> Option<(f64, f64, f64, f64, f64, f64)> {
    match comp_idx {
        IDX_CO2 => Some((304.2, 0.0408481264767325, -0.07898609940571702, 0.050275716341630625,
                         0.006766107309550113, -0.007355686498441823)),
        IDX_H2S => Some((373.2, 0.04743797774228006, -0.09161393978798885, 0.049823955842864476, 0.0, 0.0)),
        IDX_CH4 => Some((190.6, 0.1343908999845734, -0.13425828898747677, 0.040842761446578946, 0.0, 0.0)),
        IDX_N2 => Some((126.1, 0.2194835079768253, -0.14859950898753865, 0.030616725180758284, 0.0, 0.0)),
        IDX_H2 => Some((33.145, 0.37671225100706046, -0.06473044235781623, 0.0031554140401859033, 0.0, 0.0)),
        IDX_HE => Some((5.1953, 0.36301971362904206, -0.009639364760477095, 7.093021514306463e-05, 0.0, 0.0)),
        IDX_C2H6 => Some((305.4, 0.07776955167556454, -0.12336029130715784, 0.06258661601243726, 0.0, 0.0)),
        IDX_C3H8 => Some((369.8, 0.05805669163435654, -0.11105470117150809, 0.0692266175055278, 0.0, 0.0)),
        IDX_NC4H10 => Some((425.2, 0.04564160593297639, -0.10008931089682649, 0.07292670218334689, 0.0, 0.0)),
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
    0.468,  // He: H2 value as a surrogate (no He water-content data); unvalidated
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
#[allow(dead_code)]  // Sechenov salting path, not currently reached
pub fn sw_equation_8_ks(t_c: f64, tb_k: f64) -> f64 {
    let t_f = t_c * 9.0 / 5.0 + 32.0;
    0.13163 + 4.45e-4 * tb_k - 7.692e-4 * t_f
        + 2.6614e-6 * t_f * t_f
        - 2.612e-9 * t_f * t_f * t_f
}

/// Get Sechenov coefficient ks for a given gas.
///
/// For CO2 and H2S, specialized models are used on the Python side. This Rust implementation provides the S&W Eq 8 fallback
/// plus the N2 +0.02 offset.
#[allow(dead_code)]  // Sechenov salting path, not currently reached
pub fn get_sechenov_ks_sw(comp_idx: usize, t_k: f64) -> f64 {
    if comp_idx >= NUM_COMPONENTS || comp_idx == IDX_H2O {
        return 0.0;
    }
    let t_c = t_k - 273.15;
    let tb_k = COMPONENT_DB[comp_idx].tb;
    let mut ks = sw_equation_8_ks(t_c, tb_k);
    // N2: empirical +0.02 offset
    if comp_idx == IDX_N2 {
        ks += 0.02;
    }
    ks
}

// =============================================================================
// Full kij matrix construction
// =============================================================================

/// Build the N x N kij matrix for one phase of the procedure-1 flash.
///
/// `aqueous`: true for the aqueous (liquid-root) phase, where gas-water pairs
/// take the freshwater kij_AQ plus the embedded salinity delta_kij; false for
/// the non-aqueous (vapour-root) phase, where they take kij_NA. Gas-gas pairs
/// are the same in both. Returns a flattened row-major N x N matrix.
pub fn build_kij_matrix(
    comp_indices: &[usize],
    t_k: f64,
    aqueous: bool,
    salinity_molal: f64,
) -> Vec<f64> {
    let nc = comp_indices.len();
    let mut kij = vec![0.0; nc * nc];

    for i in 0..nc {
        for j in (i + 1)..nc {
            let idx_i = comp_indices[i];
            let idx_j = comp_indices[j];

            let val = if idx_i == IDX_H2O || idx_j == IDX_H2O {
                let gas_idx = if idx_i == IDX_H2O { idx_j } else { idx_i };
                if aqueous {
                    get_kij_aq_default(gas_idx, t_k)
                        + calc_embedded_delta_kij(gas_idx, t_k, salinity_molal)
                } else {
                    get_kij_na(gas_idx)
                }
            } else {
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
        // Procedure-1 CH4 at 373.15 K: Tr = 1.957765, kij = -0.10044...
        let tr = 373.15 / 190.60;
        let expect = (-2.1976782374118784 + tr) / (0.7150153669761119 + 0.8945029639230546 * tr);
        assert!((get_kij_aq_default(IDX_CH4, 373.15) - expect).abs() < 1e-15);
        assert!(get_kij_aq_default(IDX_CH4, 373.15) < 0.0);
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
    fn test_kij_he() {
        let tr = 323.15 / 5.1953;
        let expect = (-71.61550493445799 + tr) / (-0.1644070191372901 + 0.5017371389503704 * tr);
        assert!((get_kij_aq_default(IDX_HE, 323.15) - expect).abs() < 1e-15);
        assert!((get_kij_na(IDX_HE) - 0.468).abs() < 1e-15);
        assert!(calc_embedded_delta_kij(IDX_HE, 323.15, 1.0) > 0.0);
    }

    #[test]
    fn test_kij_matrix_phases() {
        // Aqueous and non-aqueous matrices differ only in gas-water pairs.
        let comp = [IDX_H2O, IDX_CH4, IDX_CO2];
        let aq = build_kij_matrix(&comp, 350.0, true, 1.0);
        let na = build_kij_matrix(&comp, 350.0, false, 1.0);
        let expect_aq = get_kij_aq_default(IDX_CO2, 350.0) + calc_embedded_delta_kij(IDX_CO2, 350.0, 1.0);
        assert!((aq[2] - expect_aq).abs() < 1e-15);
        assert!((na[2] - get_kij_na(IDX_CO2)).abs() < 1e-15);
        assert_eq!(aq[5], na[5]); // CH4-CO2
        assert_eq!(aq[5], get_gas_gas_bip(IDX_CH4, IDX_CO2));
    }

    #[test]
    fn test_sw_eq8_ks() {
        // Test at 25C for H2 (Tb=20.3 K)
        let ks = sw_equation_8_ks(25.0, 20.3);
        assert!(ks > 0.0);
        assert!(ks < 0.5);
    }
}
