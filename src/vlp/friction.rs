/// Serghides explicit Fanning friction factor.
/// Direct port of nodal.py lines 635-647.

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
    if re <= 2100.0 {
        return 16.0 / re;
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
