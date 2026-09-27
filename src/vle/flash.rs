//! Flash calculation engine: flash_tp.
//!
//! Implements the Curtis Whitson (Feb 2026) dual-flash scheme:
//!   Flash 1: All gas-water BIPs = kij_AQ → take AQUEOUS phase (gas solubilities)
//!   Flash 2: All gas-water BIPs = kij_NA → take NON-AQUEOUS phase (water content)
//!   True K-values: K_i = y_i(Flash 2) / x_i(Flash 1)

use crate::vle::alpha::{alpha_standard_pr, alpha_water_mc3, alpha_water_soreide};
use crate::vle::bip::{build_kij_matrix, Framework};
use crate::vle::components::*;
use crate::vle::fugacity::calc_fugacity_fast;
use crate::vle::k_init::sw_kvalue_init;
use crate::vle::rachford_rice::solve_rachford_rice;

/// Trivial-solution test for `flash_tp`: max|ln K| below this means both
/// labelled phases collapsed onto one composition (K -> 1). Gas-water K-values
/// sit orders of magnitude from 1, so a genuine split never comes near it.
/// Same value as Python `_lib_vle_engine.TRIVIAL_LNK_TOL`.
const TRIVIAL_LNK_TOL: f64 = 1e-6;

/// Tangent-plane stability (Michelsen 1982, Fluid Phase Equilib. 9:1-19): a
/// non-trivial stationary point with sum(Y) above 1 + this makes the feed
/// unstable. Same value as Python `_lib_vle_engine.STABILITY_TM_TOL`.
const STABILITY_TM_TOL: f64 = 1e-8;

/// Precomputed EOS quantities that are constant across SS iterations (T,P fixed).
struct EosPrecomputed {
    ai_dim: Vec<f64>,  // Dimensionless Ai = ai*P/(RT)^2
    bi_dim: Vec<f64>,  // Dimensionless Bi = bi*P/(RT)
    sqrt_ai: Vec<f64>, // sqrt(Ai)
    onemk: Vec<f64>,   // Flattened (1 - kij) matrix, nc x nc
}

/// Calculate alpha for all components. The water alpha follows the framework:
/// Soreide-Whitson (salinity-dependent) for the default framework, and
/// Mathias-Copeman 3-parameter for mc3.
fn calc_alpha(
    comp_indices: &[usize],
    t_k: f64,
    tc: &[f64],
    framework: Framework,
    salinity_molal: f64,
) -> Vec<f64> {
    let nc = comp_indices.len();
    let mut alpha = vec![0.0; nc];
    for i in 0..nc {
        if comp_indices[i] == IDX_H2O {
            let tr_w = t_k / tc[i];
            alpha[i] = match framework {
                Framework::Mc3 => alpha_water_mc3(tr_w),
                Framework::Default => alpha_water_soreide(tr_w, salinity_molal),
            };
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
    framework: Framework,
    salinity_molal: f64,
) -> EosPrecomputed {
    let nc = comp_indices.len();
    let alpha = calc_alpha(comp_indices, t_k, tc, framework, salinity_molal);
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
/// (V, x, y, converged) - vapor fraction, liquid comp, vapor comp, convergence flag.
/// When successive substitution collapses onto the trivial solution (all
/// K -> 1), a tangent-plane stability test decides: a stable feed returns
/// single phase, converged (V = 0 with x = z when water is the majority
/// component, else V = 1 with y = z; the other composition is the incipient
/// stationary point, or the feed when none exists); an unstable feed is
/// re-flashed from the trial's K. `converged` is false only if the result is
/// still trivial or the test is inconclusive. Port of Python `flash_tp`.
pub fn flash_tp(
    t_k: f64,
    p_pa: f64,
    z: &[f64],
    comp_indices: &[usize],
    mode_aq: bool,
    gamma: &[f64],
    framework: Framework,
    salinity_molal: f64,
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
    let kij_flat = build_kij_matrix(comp_indices, t_k, mode_aq, framework, salinity_molal);

    // Build component property arrays
    let tc: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].tc).collect();
    let pc: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].pc).collect();
    let omega: Vec<f64> = comp_indices.iter().map(|&i| COMPONENT_DB[i].omega).collect();

    // Precompute EOS quantities
    let eos = precompute_eos(
        comp_indices, &tc, &pc, t_k, p_pa, &kij_flat, framework, salinity_molal,
    );

    // Find water index in the component array
    let iw = comp_indices.iter().position(|&i| i == IDX_H2O).unwrap_or(0);

    // Initialize K-values (SW-specific + gamma for initial estimate)
    let k_init = sw_kvalue_init(comp_indices, &tc, &pc, &omega, t_k, p_pa);
    let mut k0: Vec<f64> = k_init
        .iter()
        .zip(gamma.iter())
        .map(|(&ki, &gi)| ki * gi)
        .collect();

    // Water K should be << 1 in gas-water systems
    k0[iw] = k0[iw].min(0.01);

    let (mut k_vals, mut converged) = ss_iterate(&z_norm, &k0, &eos, gamma, nc, max_iter, tol);

    // SS collapsed onto K = 1: both labelled phases took one composition, so
    // V says nothing. Settle it with a stability test instead.
    if is_trivial_k(&k_vals) {
        match resolve_trivial(&z_norm, &k0, &eos, gamma, iw, max_iter, tol) {
            Some(Verdict::Stable { feed_liquid, incipient }) => {
                return if feed_liquid {
                    (0.0, z_norm, incipient, true)
                } else {
                    (1.0, incipient, z_norm, true)
                };
            }
            Some(Verdict::Split(k_seed)) => {
                (k_vals, converged) = ss_iterate(&z_norm, &k_seed, &eos, gamma, nc, max_iter, tol);
            }
            None => {}
        }
    }

    // Final compositions with converged K
    let (v, mut x, mut y) = solve_rachford_rice(&z_norm, &k_vals);
    clip_and_normalize(&mut x);
    clip_and_normalize(&mut y);

    let converged = converged && !is_trivial_k(&k_vals);

    (v, x, y, converged)
}

/// True when every K is within `TRIVIAL_LNK_TOL` of 1 in log terms.
fn is_trivial_k(k: &[f64]) -> bool {
    k.iter().map(|k| k.ln().abs()).fold(0.0_f64, f64::max) < TRIVIAL_LNK_TOL
}

/// Damped successive substitution from `k_start`. Returns (K, converged).
fn ss_iterate(
    z: &[f64],
    k_start: &[f64],
    eos: &EosPrecomputed,
    gamma: &[f64],
    nc: usize,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, bool) {
    let mut k_vals = k_start.to_vec();
    let mut converged = false;
    let mut damp: f64 = 0.5;
    let mut prev_err: f64 = f64::INFINITY;

    for _it in 0..max_iter {
        // Robust RR solver (Nielsen & Lia 2022)
        let (_v, mut x, mut y) = solve_rachford_rice(z, &k_vals);

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
        let err = k_vals
            .iter()
            .zip(k_new.iter())
            .map(|(&ko, &kn)| (kn / ko - 1.0).abs())
            .fold(0.0_f64, f64::max);

        if err < tol {
            converged = true;
            k_vals = k_new;
            break;
        }

        // Adaptive damping: accelerate when converging, brake when stalling
        if err < prev_err {
            damp = (damp * 1.1).min(0.95);
        } else {
            damp = (damp * 0.5).max(0.1);
        }
        prev_err = err;

        for i in 0..nc {
            k_vals[i] *= (k_new[i] / k_vals[i]).powf(damp);
        }
    }
    (k_vals, converged)
}

/// Outcome of the stability test after a trivial SS exit.
enum Verdict {
    /// Feed is one phase; `incipient` is the opposite-phase stationary point
    /// (the feed itself when that trial collapsed).
    Stable { feed_liquid: bool, incipient: Vec<f64> },
    /// Feed splits; K seeded from the unstable trial.
    Split(Vec<f64>),
}

/// Tangent-plane stationary point by direct substitution (Michelsen 1982):
/// ln Y_i = ln f_i(z) - ln phi_i(y), less ln gamma_i for a liquid trial.
/// Returns (sum Y, normalised y), or None when phi is non-finite or the
/// search does not settle. Port of Python `_stability_trial`.
#[allow(clippy::too_many_arguments)]
fn stability_trial(
    lnf_feed: &[f64],
    seed: &[f64],
    liquid: bool,
    eos: &EosPrecomputed,
    ln_gamma: &[f64],
    nc: usize,
    max_iter: usize,
    tol: f64,
) -> Option<(f64, Vec<f64>)> {
    let seed_sum: f64 = seed.iter().sum();
    let mut ln_y: Vec<f64> = seed.iter().map(|&s| (s / seed_sum).max(1e-300).ln()).collect();
    for _ in 0..max_iter {
        let big_y: Vec<f64> = ln_y.iter().map(|v| v.exp()).collect();
        let y_sum: f64 = big_y.iter().sum();
        let y: Vec<f64> = big_y.iter().map(|v| v / y_sum).collect();
        let phi = calc_fugacity_fast(&y, &eos.ai_dim, &eos.bi_dim, &eos.sqrt_ai, &eos.onemk, nc, liquid);
        if !phi.iter().all(|p| p.is_finite()) {
            return None;
        }
        let mut step = 0.0_f64;
        for i in 0..nc {
            let mut v = lnf_feed[i] - phi[i].ln();
            if liquid {
                v -= ln_gamma[i];
            }
            step = step.max((v - ln_y[i]).abs());
            ln_y[i] = v;
        }
        if step < tol {
            let big_y: Vec<f64> = ln_y.iter().map(|v| v.exp()).collect();
            let y_sum: f64 = big_y.iter().sum();
            return Some((y_sum, big_y.iter().map(|v| v / y_sum).collect()));
        }
    }
    None
}

/// Stability test on the feed after SS collapsed onto K = 1. The feed is
/// taken as aqueous (liquid root, with gamma) when water is the majority
/// component, else as the non-aqueous phase (vapor root). Vapor-like trials
/// seeded from K0*z and from the gas components, and liquid-like trials seeded
/// from z/K0 and from pure water, search for a second phase; any unstable
/// trial wins. Port of Python `_resolve_trivial`.
fn resolve_trivial(
    z: &[f64],
    k0: &[f64],
    eos: &EosPrecomputed,
    gamma: &[f64],
    iw: usize,
    max_iter: usize,
    tol: f64,
) -> Option<Verdict> {
    let nc = z.len();
    let feed_liquid = z[iw] > 0.5;
    let ln_gamma: Vec<f64> = gamma.iter().map(|g| g.ln()).collect();
    let zc: Vec<f64> = z.iter().map(|&v| v.max(1e-300)).collect();
    let phi_z = calc_fugacity_fast(&zc, &eos.ai_dim, &eos.bi_dim, &eos.sqrt_ai, &eos.onemk, nc, feed_liquid);
    if !phi_z.iter().all(|p| p.is_finite()) {
        return None;
    }
    let lnf: Vec<f64> = (0..nc)
        .map(|i| zc[i].ln() + phi_z[i].ln() + if feed_liquid { ln_gamma[i] } else { 0.0 })
        .collect();

    // Wilson-type seeds alone can start water-dominated and collapse onto the
    // feed even when a split exists, so each trial phase also gets a near-pure
    // seed: gas-only for the vapor, water-only for the liquid.
    let gas_seed: Vec<f64> = (0..nc).map(|i| if i == iw { 0.0 } else { zc[i] } + 1e-10).collect();
    let water_seed: Vec<f64> = (0..nc).map(|i| if i == iw { 1.0 } else { 1e-10 }).collect();
    let wilson_vap: Vec<f64> = (0..nc).map(|i| k0[i] * zc[i]).collect();
    let wilson_liq: Vec<f64> = (0..nc).map(|i| zc[i] / k0[i]).collect();
    let seeds: [(bool, &[f64]); 4] = [
        (false, &wilson_vap),
        (false, &gas_seed),
        (true, &wilson_liq),
        (true, &water_seed),
    ];

    let mut incipient = z.to_vec();
    for (liquid, seed) in seeds {
        let (sum_y, y) = stability_trial(&lnf, seed, liquid, eos, &ln_gamma, nc, max_iter, tol)?;
        let k_trial: Vec<f64> = (0..nc).map(|i| y[i] / zc[i]).collect();
        if is_trivial_k(&k_trial) {
            continue;
        }
        if sum_y > 1.0 + STABILITY_TM_TOL {
            let k_seed = if liquid { k_trial.iter().map(|k| 1.0 / k).collect() } else { k_trial };
            return Some(Verdict::Split(k_seed));
        }
        if liquid != feed_liquid {
            incipient = y;
        }
    }
    Some(Verdict::Stable { feed_liquid, incipient })
}

/// Clip compositions to [1e-15, inf) and normalize.
fn clip_and_normalize(comp: &mut [f64]) {
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
            Framework::Mc3,
            0.0,    // freshwater
            200,
            1e-10,
        );

        assert!(converged, "Flash should converge");
        assert!((0.0..=1.0).contains(&v), "V should be in [0,1], got {}", v);
        assert_eq!(x.len(), 2);
        assert_eq!(y.len(), 2);

        // Water should dominate liquid phase
        assert!(x[0] > 0.9, "Liquid should be mostly water");
        // Gas should dominate vapor phase
        assert!(y[1] > 0.5, "Vapor should be mostly gas");
    }
}
