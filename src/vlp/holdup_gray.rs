/// Gray holdup and effective roughness.
/// Direct port of nodal.py lines 1304-1349.

use super::pvt_helpers::clamp;

const G_FT: f64 = 32.174;

/// Gray liquid holdup.
pub fn gray_liquid_holdup(
    v_m: f64, rho_l: f64, rho_g: f64, sigma: f64,
    diam: f64, lambda_l: f64, _p_sys: f64,
) -> f64 {
    if lambda_l <= 0.0 {
        return 0.0;
    }
    if lambda_l >= 1.0 {
        return 1.0;
    }
    if v_m < 1e-10 {
        return lambda_l;
    }

    let sigma_lbf_ft = sigma * 6.852e-5;
    if sigma_lbf_ft <= 0.0 {
        return lambda_l;
    }

    let rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l);
    if rho_ns <= 0.0 {
        return lambda_l;
    }

    let rv = lambda_l;
    let drho = (rho_l - rho_g).max(0.1);
    let n1 = v_m * v_m * rho_ns / (G_FT * diam * drho);
    let n2 = G_FT * diam * diam * drho / sigma_lbf_ft;

    let a = 0.0814 * (1.0 - 0.0554 * (1.0 + 730.0 * rv / (1.0 + rv)).ln());
    let b = 0.0554;

    if n1 <= 1e-15 {
        return lambda_l;
    }

    let argument = (a + b * n2) / n1;
    let f_e = if argument > 0.0 {
        -0.0554 * argument.ln()
    } else {
        0.0
    };

    let hl = 1.0 - (1.0 - lambda_l) * f_e.exp();
    clamp(hl, lambda_l, 1.0)
}

/// Gray effective roughness (ft).
pub fn gray_effective_roughness(
    rough_dry: f64, sigma: f64, rho_ns: f64, v_m: f64, lambda_l: f64,
) -> f64 {
    let sigma_lbf_ft = sigma * 6.852e-5;
    if v_m < 1e-10 || rho_ns <= 0.0 || sigma_lbf_ft <= 0.0 {
        return rough_dry;
    }
    let ke = rho_ns * v_m * v_m;
    let r1 = 28.5 * sigma_lbf_ft / ke.max(1e-10);
    let r2 = ke * lambda_l / sigma_lbf_ft;
    let film = r1 * (1.0 - (-r2).exp());
    rough_dry + film.max(0.0)
}
