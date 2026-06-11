/// Gray (1974) holdup and effective roughness, API 14B form.
/// Direct port of nodal.py _gray_liquid_holdup and _gray_effective_roughness.

use super::constants::*;
use super::pvt_helpers::clamp;

/// Gray (1974) liquid holdup, API 14B form.
///
/// R = vsl/vsg, fg = (1 - exp(A)) / (R + 1), HL = 1 - fg with
/// A = -2.314 * (NV * (1 + 205/ND))^B and
/// B = 0.0814 * (1 - 0.0554 * ln(1 + 730*R/(R+1))).
/// NV and ND use sigma in dyne/cm; 453.592 converts dyne/cm to lbm/s^2.
pub fn gray_liquid_holdup(
    v_sl: f64, v_sg: f64, rho_l: f64, rho_g: f64,
    sigma: f64, diam: f64, lambda_l: f64,
) -> f64 {
    if lambda_l <= 0.0 {
        return 0.0;
    }
    if lambda_l >= 1.0 {
        return 1.0;
    }
    if v_sg < 1e-10 {
        return 1.0;
    }
    if sigma <= 0.0 {
        return lambda_l;
    }

    let rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l);
    if rho_ns <= 0.0 {
        return lambda_l;
    }

    let v_m = v_sl + v_sg;
    let r = v_sl / v_sg;
    let drho = (rho_l - rho_g).max(0.1);
    let nv = DYNCM_PER_LBM_S2 * rho_ns.powi(2) * v_m.powi(4) / (G_FT * sigma * drho);
    let nd = DYNCM_PER_LBM_S2 * G_FT * drho * diam.powi(2) / sigma;

    let b = GRAY_B1 * (1.0 - GRAY_B2 * (1.0 + GRAY_B3 * r / (r + 1.0)).ln());
    let a = -GRAY_A_COEF * (nv * (1.0 + GRAY_ND_COEF / nd)).powf(b);
    let fg = (1.0 - a.exp()) / (r + 1.0);
    clamp(1.0 - fg, lambda_l, 1.0)
}

/// Gray (1974) effective (wet film) roughness in ft, API 14B form.
///
/// ke0 = 28.5 * sigma / (453.592 * rho_ns * vm^2) with sigma in dyne/cm.
/// For R = vsl/vsg >= 0.007 use ke0; below, interpolate from the dry pipe
/// roughness. Result floored at 2.77e-5 ft.
pub fn gray_effective_roughness(
    rough_dry: f64, sigma: f64, rho_ns: f64, v_sl: f64, v_sg: f64,
) -> f64 {
    let v_m = v_sl + v_sg;
    if v_m < 1e-10 || rho_ns <= 0.0 || sigma <= 0.0 || v_sg < 1e-10 {
        return rough_dry;
    }
    let r = v_sl / v_sg;
    let ke0 = GRAY_ROUGH_K * sigma / (DYNCM_PER_LBM_S2 * rho_ns * v_m * v_m);
    let ke = if r >= GRAY_R_THRESH {
        ke0
    } else {
        rough_dry + r * (ke0 - rough_dry) / GRAY_R_THRESH
    };
    ke.max(GRAY_ROUGH_FLOOR)
}
