/// Hagedorn-Brown holdup correlation helpers.
/// Direct port of nodal.py HB polynomial holdup (CNL, YLONSI, SI, Orkiszewski).

use super::pvt_helpers::clamp;

const LN10: f64 = std::f64::consts::LN_10;

#[inline]
fn log10_safe(x: f64) -> f64 {
    if x <= 0.0 { -30.0 } else { x.ln() / LN10 }
}

/// Compute HB liquid holdup from dimensionless numbers.
/// Returns (yl, nvl, nvg, nd).
pub fn hb_holdup(
    ul: f64, ugas: f64, tid: f64, lsg: f64, ift: f64,
    mul: f64, p_avg: f64,
) -> (f64, f64, f64, f64) {
    let nvl = 1.938 * ul * (62.4 * lsg / ift).powf(0.25);
    let nvg = 1.938 * ugas * (62.4 * lsg / ift).powf(0.25);
    let nd = 120.872 * tid / 12.0 * (62.4 * lsg / ift).powf(0.5);

    let nl = if mul > 0.0 {
        let v = 0.15726 * mul * (1.0 / (62.4 * lsg * ift.powi(3))).powf(0.25);
        if v <= 0.0 { 1e-8 } else { v }
    } else {
        1e-8
    };

    let x1 = log10_safe(nl) + 3.0;
    let cnl = 10.0_f64.powf(
        -2.69851 + 0.1584095 * x1 - 0.5509976 * x1 * x1
            + 0.5478492 * x1.powi(3) - 0.1219458 * x1.powi(4),
    );

    let nvg_safe = nvg.max(1e-10);
    let f1 = nvl * p_avg.powf(0.1) * cnl / (nvg_safe.powf(0.575) * 14.7_f64.powf(0.1) * nd);

    let lf1 = log10_safe(f1) + 6.0;
    let mut ylonsi = -0.10306578 + 0.617774 * lf1 - 0.632946 * lf1.powi(2)
        + 0.29598 * lf1.powi(3) - 0.0401 * lf1.powi(4);
    ylonsi = ylonsi.max(0.0);

    let mut f2 = nvg * nl.powf(0.38) / nd.powf(2.14);
    let idex: f64 = if f2 >= 0.012 { 1.0 } else { -1.0 };
    f2 = (1.0 - idex) / 2.0 * 0.012 + (1.0 + idex) / 2.0 * f2;

    let si = if f2 <= 0.001 {
        1.0
    } else {
        0.9116257 - 4.821756 * f2 + 1232.25 * f2 * f2
            - 22253.58 * f2.powi(3) + 116174.3 * f2.powi(4)
    };

    let yl = clamp(si * ylonsi, 0.0, 1.0);

    (yl, nvl, nvg, nd)
}

/// Orkiszewski bubble flow correction applied to yl.
pub fn orkiszewski_correction(yl: f64, ugas: f64, ul: f64, tid: f64) -> f64 {
    let vm = ugas + ul;
    let vs = 0.8_f64;
    let d_ft = tid / 12.0;
    let mut lb = 1.071 - 0.2218 * vm * vm / d_ft;
    if lb < 0.13 {
        lb = 0.13;
    }
    if vm > 0.0 {
        let b_ratio = ugas / vm;
        if b_ratio < lb {
            let disc = (1.0 + vm / vs).powi(2) - 4.0 * ugas / vs;
            if disc >= 0.0 {
                let yl_new = 1.0 - 0.5 * (1.0 + vm / vs - disc.sqrt());
                return clamp(yl_new, 0.0, 1.0);
            }
        }
    }
    yl
}
