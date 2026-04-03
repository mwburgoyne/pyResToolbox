/// Flash calculation engine: flash_tp and calc_equilibrium.
///
/// Implements the Curtis Whitson (Feb 2026) dual-flash scheme:
///   Flash 1: All gas-water BIPs = kij_AQ → take AQUEOUS phase (gas solubilities)
///   Flash 2: All gas-water BIPs = kij_NA → take NON-AQUEOUS phase (water content)
///   True K-values: K_i = y_i(Flash 2) / x_i(Flash 1)

use crate::vle::alpha::{alpha_standard_pr, alpha_water_mc3};
use crate::vle::bip::{build_kij_matrix, get_sechenov_ks_sw};
use crate::vle::components::*;
use crate::vle::fugacity::calc_fugacity_fast;
use crate::vle::k_init::sw_kvalue_init;
use crate::vle::rachford_rice::solve_rachford_rice;

/// Precomputed EOS quantities that are constant across SS iterations (T,P fixed).
struct EosPrecomputed {
    ai_dim: Vec<f64>,  // Dimensionless Ai = ai*P/(RT)^2
    bi_dim: Vec<f64>,  // Dimensionless Bi = bi*P/(RT)
    sqrt_ai: Vec<f64>, // sqrt(Ai)
    onemk: Vec<f64>,   // Flattened (1 - kij) matrix, nc x nc
}

/// Calculate alpha for all components using proposed framework (MC-3 for water).
fn calc_alpha_proposed(comp_indices: &[usize], t_k: f64, tc: &[f64]) -> Vec<f64> {
    let nc = comp_indices.len();
    let mut alpha = vec![0.0; nc];
    for i in 0..nc {
        if comp_indices[i] == IDX_H2O {
            let tr_w = t_k / tc[i];
            alpha[i] = alpha_water_mc3(tr_w);
        } else {
            let tr = t_k / tc[i];
            alpha[i] = alpha_standard_pr(tr, COMPONENT_DB[comp_indices[i]].omega);
        }
    }
    alpha
}

/// Calculate ai, bi from alpha, Tc, Pc.
fn calc_ai_bi(alpha: &[f64], tc: &[f64], pc: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let nc = alpha.len();
    let mut ai = vec![0.0; nc];
    let mut bi = vec![0.0; nc];
    for i in 0..nc {
        ai[i] = OMEGA_A * (R_GAS * tc[i]).powi(2) * alpha[i] / pc[i];
        bi[i] = OMEGA_B * R_GAS * tc[i] / pc[i];
    }
    (ai, bi)
}

/// Precompute EOS quantities for given T, P, and kij matrix.
fn precompute_eos(
    comp_indices: &[usize],
    tc: &[f64],
    pc: &[f64],
    t_k: f64,
    p_pa: f64,
    kij_flat: &[f64],
) -> EosPrecomputed {
    let nc = comp_indices.len();
    let alpha = calc_alpha_proposed(comp_indices, t_k, tc);
    let (ai, bi) = calc_ai_bi(&alpha, tc, pc);

    let rt = R_GAS * t_k;
    let rt2 = rt * rt;
    let mut ai_dim = vec![0.0; nc];
    let mut bi_dim = vec![0.0; nc];
    let mut sqrt_ai_dim = vec![0.0; nc];

    for i in 0..nc {
        ai_dim[i] = ai[i] * p_pa / rt2;
        bi_dim[i] = bi[i] * p_pa / rt;
        sqrt_ai_dim[i] = ai_dim[i].sqrt();
    }

    // onemk = 1.0 - kij_matrix (flattened)
    let onemk: Vec<f64> = kij_flat.iter().map(|&k| 1.0 - k).collect();

    EosPrecomputed {
        ai_dim,
        bi_dim,
        sqrt_ai: sqrt_ai_dim,
        onemk,
    }
}

/// Full TP flash using successive substitution with robust RR solver.
///
/// # Arguments
/// * `t_k` - Temperature in Kelvin
/// * `p_pa` - Pressure in Pascal
/// * `z` - Feed composition (will be normalized)
/// * `comp_indices` - Component index array
/// * `mode_aq` - true for AQ mode, false for NA mode
/// * `gamma` - Activity coefficient array (length nc). Use all-ones for freshwater.
/// * `max_iter` - Maximum SS iterations
/// * `tol` - Convergence tolerance on K-values
///
/// # Returns
/// (V, x, y, converged) - vapor fraction, liquid comp, vapor comp, convergence flag
pub fn flash_tp(
    t_k: f64,
    p_pa: f64,
    z: &[f64],
    comp_indices: &[usize],
    mode_aq: bool,
    gamma: &[f64],
    max_iter: usize,
    tol: f64,
) -> (f64, Vec<f64>, Vec<f64>, bool) {
    let nc = comp_indices.len();
    assert_eq!(z.len(), nc);
    assert_eq!(gamma.len(), nc);

    // Normalize z
    let z_sum: f64 = z.iter().sum();
    let z_norm: Vec<f64> = z.iter().map(|&x| x / z_sum).collect();

    // Build kij matrix
    let kij_flat = build_kij_matrix(comp_indices, t_k, mode_aq);

    // Build component property arrays
    let tc: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].tc).collect();
    let pc: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].pc).collect();
    let omega: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].omega).collect();

    // Precompute EOS quantities
    let eos = precompute_eos(comp_indices, &tc, &pc, t_k, p_pa, &kij_flat);

    // Find water index in the component array
    let iw = comp_indices.iter().position(|&i| i == IDX_H2O).unwrap_or(0);

    // Initialize K-values (SW-specific + gamma for initial estimate)
    let k_init = sw_kvalue_init(comp_indices, &tc, &pc, &omega, t_k, p_pa);
    let mut k_vals: Vec<f64> = k_init
        .iter()
        .zip(gamma.iter())
        .map(|(&ki, &gi)| ki * gi)
        .collect();

    // Water K should be << 1 in gas-water systems
    k_vals[iw] = k_vals[iw].min(0.01);

    let mut converged = false;

    for it in 0..max_iter {
        // Robust RR solver (Nielsen & Lia 2022)
        let (_v, mut x, mut y) = solve_rachford_rice(&z_norm, &k_vals);

        // Clip and normalize
        clip_and_normalize(&mut x);
        clip_and_normalize(&mut y);

        // Fugacity coefficients
        let phi_l = calc_fugacity_fast(&x, &eos.ai_dim, &eos.bi_dim, &eos.sqrt_ai, &eos.onemk, nc, true);
        let phi_v = calc_fugacity_fast(&y, &eos.ai_dim, &eos.bi_dim, &eos.sqrt_ai, &eos.onemk, nc, false);

        // Gamma-phi K-value: K_i = gamma_i * phi_L_i / phi_V_i
        let mut k_new: Vec<f64> = Vec::with_capacity(nc);
        for i in 0..nc {
            let ki = (gamma[i] * phi_l[i] / (phi_v[i] + 1e-30)).clamp(1e-10, 1e10);
            k_new.push(ki);
        }

        // Check convergence
        let max_rel_change = k_vals
            .iter()
            .zip(k_new.iter())
            .map(|(&ko, &kn)| (kn / ko - 1.0).abs())
            .fold(0.0_f64, f64::max);

        if max_rel_change < tol {
            converged = true;
            k_vals = k_new;
            break;
        }

        // Damped successive substitution
        let damp = if it < 20 { 0.7 } else { 0.9 };
        for i in 0..nc {
            k_vals[i] *= (k_new[i] / k_vals[i]).powf(damp);
        }
    }

    // Final compositions with converged K
    let (v, mut x, mut y) = solve_rachford_rice(&z_norm, &k_vals);
    clip_and_normalize(&mut x);
    clip_and_normalize(&mut y);

    (v, x, y, converged)
}

/// Clip compositions to [1e-15, inf) and normalize.
fn clip_and_normalize(comp: &mut Vec<f64>) {
    for v in comp.iter_mut() {
        if *v < 1e-15 {
            *v = 1e-15;
        }
    }
    let sum: f64 = comp.iter().sum();
    for v in comp.iter_mut() {
        *v /= sum;
    }
}

/// Calculate Sechenov activity coefficients for gamma-phi flash.
///
/// gamma_i = 10^(ks_i * m) for gases, gamma_H2O = 1.0
///
/// Uses S&W Equation 8 for all gases (the proposed framework routes CO2/H2S
/// to specialized models on the Python side; for Rust we use the S&W Eq 8 fallback).
pub fn calc_gamma(
    comp_indices: &[usize],
    t_k: f64,
    salinity_molal: f64,
) -> Vec<f64> {
    let nc = comp_indices.len();
    let mut gamma = vec![1.0; nc];

    if salinity_molal <= 0.0 {
        return gamma;
    }

    for i in 0..nc {
        if comp_indices[i] != IDX_H2O {
            let ks = get_sechenov_ks_sw(comp_indices[i], t_k);
            gamma[i] = 10.0_f64.powf(ks * salinity_molal);
        }
    }

    gamma
}

/// Full S&W equilibrium per Curtis's scheme.
///
/// Flash 1 (kij_AQ) -> aqueous phase x
/// Flash 2 (kij_NA) -> non-aqueous phase y
/// True K-values: K_i = y_na / x_aq
///
/// # Arguments
/// * `t_k` - Temperature in Kelvin
/// * `p_pa` - Pressure in Pascal
/// * `z` - Feed composition
/// * `comp_indices` - Component index array
/// * `salinity_molal` - Salinity for gamma-phi approach
///
/// # Returns
/// FlashResult with all equilibrium data.
pub struct FlashResult {
    pub x_aq: Vec<f64>,
    pub y_na: Vec<f64>,
    pub k_true: Vec<f64>,
    pub v_aq: f64,
    pub v_na: f64,
    pub converged_aq: bool,
    pub converged_na: bool,
    pub gamma: Vec<f64>,
    /// Fugacity coefficients from AQ flash (liquid phase)
    pub fug_l: Vec<f64>,
    /// Fugacity coefficients from NA flash (vapor phase)
    pub fug_v: Vec<f64>,
}

pub fn calc_equilibrium(
    t_k: f64,
    p_pa: f64,
    z: &[f64],
    comp_indices: &[usize],
    salinity_molal: f64,
) -> FlashResult {
    let nc = comp_indices.len();

    // Normalize z
    let z_sum: f64 = z.iter().sum();
    let z_norm: Vec<f64> = z.iter().map(|&x| x / z_sum).collect();

    // Calculate gamma for gamma-phi method
    let gamma = calc_gamma(comp_indices, t_k, salinity_molal);

    // Flash 1: Aqueous phase (kij_AQ) with gamma
    let (v_aq, x_aq, _y_aq, conv_aq) = flash_tp(
        t_k,
        p_pa,
        &z_norm,
        comp_indices,
        true, // AQ mode
        &gamma,
        200,
        1e-10,
    );

    // Flash 2: Non-aqueous phase (kij_NA) — no gamma (salt stays in liquid)
    let ones = vec![1.0; nc];
    let (v_na, _x_na, y_na, conv_na) = flash_tp(
        t_k,
        p_pa,
        &z_norm,
        comp_indices,
        false, // NA mode
        &ones,
        200,
        1e-10,
    );

    // True K-values: K_i = y_i(NA) / x_i(AQ)
    let k_true: Vec<f64> = x_aq
        .iter()
        .zip(y_na.iter())
        .map(|(&xi, &yi)| {
            if xi > 1e-15 {
                yi / xi
            } else {
                1e10
            }
        })
        .collect();

    // Compute fugacity coefficients for the final compositions
    let tc: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].tc).collect();
    let pc: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].pc).collect();

    // AQ flash fugacity (liquid phase)
    let kij_aq = build_kij_matrix(comp_indices, t_k, true);
    let eos_aq = precompute_eos(comp_indices, &tc, &pc, t_k, p_pa, &kij_aq);
    let fug_l = calc_fugacity_fast(
        &x_aq, &eos_aq.ai_dim, &eos_aq.bi_dim, &eos_aq.sqrt_ai, &eos_aq.onemk, nc, true,
    );

    // NA flash fugacity (vapor phase)
    let kij_na = build_kij_matrix(comp_indices, t_k, false);
    let eos_na = precompute_eos(comp_indices, &tc, &pc, t_k, p_pa, &kij_na);
    let fug_v = calc_fugacity_fast(
        &y_na, &eos_na.ai_dim, &eos_na.bi_dim, &eos_na.sqrt_ai, &eos_na.onemk, nc, false,
    );

    FlashResult {
        x_aq,
        y_na,
        k_true,
        v_aq,
        v_na,
        converged_aq: conv_aq,
        converged_na: conv_na,
        gamma,
        fug_l,
        fug_v,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flash_tp_two_component() {
        // Simple H2O + CH4 flash at moderate conditions
        let comp = [IDX_H2O, IDX_CH4];
        let z = [0.95, 0.05];
        let gamma = [1.0, 1.0];

        let (v, x, y, converged) = flash_tp(
            373.15,   // 100C
            100.0e5,  // 100 bar
            &z,
            &comp,
            true,   // AQ mode
            &gamma,
            200,
            1e-10,
        );

        assert!(converged, "Flash should converge");
        assert!(v >= 0.0 && v <= 1.0, "V should be in [0,1], got {}", v);
        assert_eq!(x.len(), 2);
        assert_eq!(y.len(), 2);

        // Water should dominate liquid phase
        assert!(x[0] > 0.9, "Liquid should be mostly water");
        // Gas should dominate vapor phase
        assert!(y[1] > 0.5, "Vapor should be mostly gas");
    }

    #[test]
    fn test_calc_equilibrium_freshwater() {
        let comp = [IDX_H2O, IDX_CH4];
        let z = [0.95, 0.05];

        let result = calc_equilibrium(373.15, 100.0e5, &z, &comp, 0.0);

        assert!(result.converged_aq, "AQ flash should converge");
        assert!(result.converged_na, "NA flash should converge");
        assert_eq!(result.x_aq.len(), 2);
        assert_eq!(result.y_na.len(), 2);
        assert_eq!(result.k_true.len(), 2);
    }

    #[test]
    fn test_calc_equilibrium_brine() {
        let comp = [IDX_H2O, IDX_CH4];
        let z = [0.95, 0.05];

        let result_fresh = calc_equilibrium(373.15, 100.0e5, &z, &comp, 0.0);
        let result_brine = calc_equilibrium(373.15, 100.0e5, &z, &comp, 1.5);

        // Brine should have LESS dissolved gas than freshwater (salting out)
        assert!(
            result_brine.x_aq[1] < result_fresh.x_aq[1],
            "Brine should have less dissolved gas"
        );
    }

    #[test]
    fn test_calc_gamma() {
        let comp = [IDX_H2O, IDX_CH4];
        let gamma = calc_gamma(&comp, 373.15, 1.5);
        assert!((gamma[0] - 1.0).abs() < 1e-10, "Water gamma should be 1.0");
        assert!(gamma[1] > 1.0, "Gas gamma should be > 1 for salting out");
    }
}
