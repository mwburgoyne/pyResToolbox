//! Serghides explicit Fanning friction factor.
//! Direct port of nodal.py _serghides_fanning.

use super::constants::{
    CHURCHILL_RE_MAX, CH_A_K, CH_A_RE, CH_A_ROUGH, CH_B_RE, CH_LAM, DR_LOG_K, DR_ROUGH_EXP,
    DR_ROUGH_K,
};

const LN10: f64 = std::f64::consts::LN_10;

#[inline]
fn log10_safe(x: f64) -> f64 {
    if x <= 0.0 { -30.0 } else { x.ln() / LN10 }
}

/// Serghides Fanning friction factor.
pub fn serghides_fanning(re: f64, eps_d: f64) -> f64 {
    if re < 1.0 {
        return 0.0;
    }
    if re < CHURCHILL_RE_MAX {
        return churchill_fanning(re, eps_d);
    }
    let a = -2.0 * log10_safe(eps_d / 3.7 + 12.0 / re);
    let b = -2.0 * log10_safe(eps_d / 3.7 + 2.51 * a / re);
    let c = -2.0 * log10_safe(eps_d / 3.7 + 2.51 * b / re);
    let diff = c - 2.0 * b + a;
    if diff.abs() < 1e-30 {
        return 0.005;
    }
    let f_darcy = (a - (b - a).powi(2) / diff).powi(-2);
    f_darcy / 4.0
}

/// Churchill (1977) Fanning friction factor, valid for all Re; used below
/// CHURCHILL_RE_MAX (Python _churchill_fanning, nodal.py).
pub fn churchill_fanning(re: f64, eps_d: f64) -> f64 {
    let a = (CH_A_K * (1.0 / ((CH_A_RE / re).powf(0.9) + CH_A_ROUGH * eps_d)).ln()).powi(16);
    let b = (CH_B_RE / re).powi(16);
    let f_darcy = CH_LAM * ((CH_LAM / re).powi(12) + (a + b).powf(-1.5)).powf(1.0 / 12.0);
    f_darcy / 4.0
}

/// Duns and Ros (1963) large-roughness Fanning friction factor, used above
/// eps/d 0.05 (Python _duns_ros_rough_fanning, nodal.py).
pub fn duns_ros_rough_fanning(eps_d: f64) -> f64 {
    1.0 / (4.0 * log10_safe(DR_LOG_K * eps_d)).powi(2) + DR_ROUGH_K * eps_d.powf(DR_ROUGH_EXP)
}
