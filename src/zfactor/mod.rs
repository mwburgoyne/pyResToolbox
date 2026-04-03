use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use crate::critical_properties;

const DEGF2R: f64 = 459.67;
const R_GAS: f64 = 10.7316;
const MW_AIR: f64 = 28.967;
const MW_CO2: f64 = 44.01;
const MW_H2S: f64 = 34.082;
const MW_N2: f64 = 28.014;
const MW_H2: f64 = 2.016;

// =========================================================================
// DAK Z-factor (Dranchuk & Abou-Kassem 1975)
// =========================================================================

/// DAK Z-factor from reduced pressure and temperature.
#[pyfunction]
pub fn dak_zfactor(pr: f64, tr: f64) -> PyResult<f64> {
    if !pr.is_finite() || !tr.is_finite() || tr == 0.0 {
        return Ok(f64::NAN);
    }
    Ok(dak_core(pr, tr))
}

fn dak_core(pr: f64, tr: f64) -> f64 {
    let a1: f64 = 0.3265;
    let a2: f64 = -1.0700;
    let a3: f64 = -0.5339;
    let a4: f64 = 0.01569;
    let a5: f64 = -0.05165;
    let a6: f64 = 0.5475;
    let a7: f64 = -0.7361;
    let a8: f64 = 0.1844;
    let a9: f64 = 0.1056;
    let a10: f64 = 0.6134;
    let a11: f64 = 0.7210;

    let tr2 = tr * tr;
    let tr3 = tr2 * tr;
    let tr4 = tr3 * tr;
    let tr5 = tr4 * tr;

    let r1 = a1 + a2 / tr + a3 / tr3 + a4 / tr4 + a5 / tr5;
    let r2 = 0.27 * pr / tr;
    let r3 = a6 + a7 / tr + a8 / tr2;
    let r4 = a9 * (a7 / tr + a8 / tr2);
    let r5 = a10 / tr3;

    let mut rhor = (0.27 * pr / tr).max(1e-10);
    let tolerance = 1e-6;

    for _ in 0..100 {
        let r2_val = rhor * rhor;
        let r5_val = rhor * rhor * rhor * rhor * rhor;
        let exp_term = (-a11 * r2_val).exp();

        let f_val = r1 * rhor - r2 / rhor + r3 * r2_val - r4 * r5_val
            + r5 * r2_val * (1.0 + a11 * r2_val) * exp_term + 1.0;

        if f_val.abs() < tolerance {
            break;
        }

        let fp_val = r1 + r2 / (rhor * rhor)
            + 2.0 * r3 * rhor
            - 5.0 * r4 * rhor * rhor * rhor * rhor
            + 2.0 * r5 * rhor * exp_term
                * ((1.0 + 2.0 * a11 * rhor * rhor * rhor) - a11 * r2_val * (1.0 + a11 * r2_val));

        let fp_safe = if fp_val.abs() < 1e-30 { 1e-30 } else { fp_val };
        let step = f_val / fp_safe;
        let rhor_new = (rhor - step).max(1e-10);

        if (rhor - rhor_new).abs() < tolerance {
            rhor = rhor_new;
            break;
        }
        rhor = rhor_new;
    }

    0.27 * pr / (rhor * tr)
}

// =========================================================================
// Hall-Yarborough Z-factor (1973)
// =========================================================================

/// Hall-Yarborough Z-factor from reduced pressure and temperature.
#[pyfunction]
pub fn hall_yarborough_zfactor(pr: f64, tr: f64) -> PyResult<f64> {
    if !pr.is_finite() || !tr.is_finite() || tr == 0.0 {
        return Ok(f64::NAN);
    }
    Ok(hy_core(pr, tr))
}

fn hy_core(pr: f64, tr: f64) -> f64 {
    let tpr_inv = 1.0 / tr;
    let t2 = tpr_inv * tpr_inv;
    let a = 0.06125 * tpr_inv * (-1.2 * (1.0 - tpr_inv).powi(2)).exp();
    let b = tpr_inv * (14.76 - 9.76 * tpr_inv + 4.58 * t2);
    let c = tpr_inv * (90.7 - 242.2 * tpr_inv + 42.4 * t2);
    let d_exp = 2.18 + 2.82 * tpr_inv;

    // Initial guess from WYW
    let z_init = wyw_core(pr, tr);
    let mut yi = (a * pr / z_init).max(1e-10);
    let mut y = 0.01_f64;

    for _ in 0..100 {
        let rel_err = (y - yi).abs() / y.abs().max(1e-10);
        if rel_err <= 0.0005 {
            break;
        }
        let yi_safe = yi.clamp(1e-10, 0.99);
        let omy = 1.0 - yi_safe;
        let omy3 = omy * omy * omy;
        let omy4 = omy3 * omy;
        let yi2 = yi_safe * yi_safe;
        let yi3 = yi2 * yi_safe;
        let yi4 = yi3 * yi_safe;

        let f_val = (yi_safe + yi2 + yi3 - yi4) / omy3
            - a * pr - b * yi2 + c * yi_safe.powf(d_exp);
        let df_val = (1.0 + 4.0 * yi_safe + 4.0 * yi2 - 4.0 * yi3 + yi4) / omy4
            - 2.0 * b * yi_safe + c * d_exp * yi_safe.powf(d_exp - 1.0);

        let df_safe = if df_val.abs() < 1e-30 { 1e-30 } else { df_val };
        y = (yi_safe - f_val / df_safe).max(1e-10);
        yi = y;
    }

    y = y.max(1e-30);
    a * pr / y
}

// =========================================================================
// WYW Z-factor (Wang, Ye & Wu 2021) - explicit, no iteration
// =========================================================================

fn wyw_core(pr: f64, tr: f64) -> f64 {
    let a: [f64; 21] = [
        0.0, 256.41675, 7.18202, -178.5725, 182.98704, -40.74427, 2.24427,
        47.44825, 5.2852, -0.14914, 271.50446, 16.2694, -121.51728,
        167.71477, -81.73093, 20.36191, -2.1177, 124.64444, -6.74331,
        0.20897, -0.00314,
    ];

    let tr2 = tr * tr;
    let tr3 = tr2 * tr;
    let tr4 = tr3 * tr;
    let tr5 = tr4 * tr;

    let pr2 = pr * pr;
    let pr3 = pr2 * pr;
    let pr4 = pr3 * pr;
    let pr5 = pr4 * pr;

    let num = a[1] + a[2] * (1.0 + a[3] * tr + a[4] * tr2 + a[5] * tr3 + a[6] * tr4) * pr
        + a[7] * pr2 + a[8] * pr3 + a[9] * pr4;
    let den = a[10] + a[11] * (1.0 + a[12] * tr + a[13] * tr2 + a[14] * tr3 + a[15] * tr4 + a[16] * tr5) * pr
        + a[17] * pr2 + a[18] * pr3 + a[19] * pr4 + a[20] * pr5;

    num / den
}

// =========================================================================
// Full pipeline functions
// =========================================================================

/// DAK full pipeline: Sutton -> Wichert-Aziz -> DAK
#[pyfunction]
pub fn dak_zfactor_full(
    p_psia: f64,
    t_degf: f64,
    sg: f64,
    co2_frac: f64,
    h2s_frac: f64,
    n2_frac: f64,
) -> PyResult<f64> {
    let (tpc, ppc) = critical_properties::sutton_wa_internal(sg, co2_frac, h2s_frac, n2_frac)
        .map_err(|e| PyValueError::new_err(e))?;

    let t_rankine = t_degf + DEGF2R;
    let pr = p_psia / ppc;
    let tr = t_rankine / tpc;
    Ok(dak_core(pr, tr))
}

/// Hall-Yarborough full pipeline: Sutton -> Wichert-Aziz -> HY
#[pyfunction]
pub fn hall_yarborough_zfactor_full(
    p_psia: f64,
    t_degf: f64,
    sg: f64,
    co2_frac: f64,
    h2s_frac: f64,
    n2_frac: f64,
) -> PyResult<f64> {
    let (tpc, ppc) = critical_properties::sutton_wa_internal(sg, co2_frac, h2s_frac, n2_frac)
        .map_err(|e| PyValueError::new_err(e))?;

    let t_rankine = t_degf + DEGF2R;
    let pr = p_psia / ppc;
    let tr = t_rankine / tpc;
    Ok(hy_core(pr, tr))
}

// =========================================================================
// BNS Z-factor (Peng-Robinson EOS, 5-component)
// =========================================================================

// EOS parameters: [CO2=0, H2S=1, N2=2, H2=3, Gas=4]
const BNS_MWS: [f64; 5] = [44.01, 34.082, 28.014, 2.016, 0.0]; // Gas MW set at runtime
const BNS_TCS: [f64; 5] = [547.416, 672.120, 227.160, 47.430, 1.0]; // Gas Tc set at runtime
const BNS_PCS: [f64; 5] = [1069.51, 1299.97, 492.84, 187.5300, 1.0]; // Gas Pc set at runtime
const BNS_ACF: [f64; 5] = [0.12253, 0.04909, 0.037, -0.21700, -0.03899];
const BNS_VSHIFT: [f64; 5] = [-0.27607, -0.22901, -0.21066, -0.36270, -0.19076];
const BNS_OMEGAA: [f64; 5] = [0.427671, 0.436725, 0.457236, 0.457236, 0.457236];
const BNS_OMEGAB: [f64; 5] = [0.0696397, 0.0724345, 0.0777961, 0.0777961, 0.0777961];

// BIP constant matrix (5x5, row-major)
const BIP_CONST: [[f64; 5]; 5] = [
    [ 0.0,       0.248638, -0.25,     -0.247153, -0.145561],
    [ 0.248638,  0.0,      -0.204414,  0.0,       0.16852],
    [-0.25,     -0.204414,  0.0,      -0.166253, -0.108],
    [-0.247153,  0.0,      -0.166253,  0.0,      -0.0620119],
    [-0.145561,  0.16852,  -0.108,    -0.0620119, 0.0],
];

// BIP slope/Tc matrix (5x5)
const BIP_SLOPE_TC: [[f64; 5]; 5] = [
    [  0.0,          -75.64467996,  63.51120432, 89.65031832, 0.0],
    [-75.64467996,    0.0,         157.55635404,  0.0,        0.0],
    [ 63.51120432,  157.55635404,   0.0,         17.90313836, 0.0],
    [ 89.65031832,    0.0,          17.90313836,  0.0,        0.0],
    [  0.0,            0.0,          0.0,          0.0,        0.0],
];

// Gas column BIP slopes
const BIP_GAS_SLOPES: [f64; 4] = [0.276572, -0.122378, 0.0605506, 0.0427873];

fn calc_bips(deg_r: f64, tpc_hc: f64) -> [[f64; 5]; 5] {
    let mut kij = [[0.0_f64; 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            let mut slope = BIP_SLOPE_TC[i][j];
            // Gas column/row adjustments
            if i == 4 && j < 4 {
                slope = BIP_GAS_SLOPES[j] * tpc_hc;
            } else if j == 4 && i < 4 {
                slope = BIP_GAS_SLOPES[i] * tpc_hc;
            }
            kij[i][j] = BIP_CONST[i][j] + slope / deg_r;
        }
    }
    kij
}

/// Cardano cubic solver for monic cubic z^3 + c2*z^2 + c1*z + c0 = 0
/// flag: 1=max root, -1=min root
fn cardano_cubic(c2: f64, c1: f64, c0: f64, flag: i32) -> f64 {
    let p = (3.0 * c1 - c2 * c2) / 3.0;
    let q = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1 + 27.0 * c0) / 27.0;
    let disc = q * q / 4.0 + p * p * p / 27.0;

    if disc < 0.0 {
        let m = 2.0 * (-p / 3.0).sqrt();
        let qpm = (3.0 * q / (p * m)).clamp(-1.0, 1.0);
        let theta1 = qpm.acos() / 3.0;
        let r0 = m * theta1.cos() - c2 / 3.0;
        let r1 = m * (theta1 + 4.0 * std::f64::consts::PI / 3.0).cos() - c2 / 3.0;
        let r2 = m * (theta1 + 2.0 * std::f64::consts::PI / 3.0).cos() - c2 / 3.0;
        if flag == 1 {
            r0.max(r1).max(r2)
        } else {
            r0.min(r1).min(r2)
        }
    } else {
        let sqrt_disc = disc.sqrt();
        let pp = -q / 2.0 + sqrt_disc;
        let qq = -q / 2.0 - sqrt_disc;
        let pp_cr = if pp >= 0.0 { pp.cbrt() } else { -(-pp).cbrt() };
        let qq_cr = if qq >= 0.0 { qq.cbrt() } else { -(-qq).cbrt() };
        pp_cr + qq_cr - c2 / 3.0
    }
}

/// Halley's method for finding max root of cubic z^3 + c2*z^2 + c1*z + c0 = 0
fn halley_cubic(c2: f64, c1: f64, c0: f64) -> f64 {
    let mut z = -c2 / 3.0;
    let f0 = z * z * z + c2 * z * z + c1 * z + c0;
    if f0 < 0.0 {
        z += 1.0;
    }

    for _ in 0..50 {
        let f = z * z * z + c2 * z * z + c1 * z + c0;
        let fp = 3.0 * z * z + 2.0 * c2 * z + c1;
        let fpp = 6.0 * z + 2.0 * c2;

        let fp_safe = if fp.abs() < 1e-30 { 1e-30 } else { fp };
        let dz = f / fp_safe;
        let denom = fp_safe - 0.5 * dz * fpp;
        let denom_safe = if denom.abs() < 1e-30 { 1e-30 } else { denom };
        let dz_final = f / denom_safe;
        z -= dz_final;
        if dz_final.abs() < 1e-12 {
            break;
        }
    }

    // Verify convergence
    let f = z * z * z + c2 * z * z + c1 * z + c0;
    if f.abs() > 1e-6 {
        // Fallback to Cardano
        return cardano_cubic(c2, c1, c0, 1);
    }
    z
}

fn bns_zfactor_core(
    p_psia: f64,
    deg_r: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
    tpc_hc: f64,
    ppc_hc: f64,
) -> f64 {
    let zi = [co2, h2s, n2, h2, 1.0 - co2 - h2s - n2 - h2];

    // Set up component properties
    let mut tcs = BNS_TCS;
    let mut pcs = BNS_PCS;
    tcs[4] = tpc_hc;
    pcs[4] = ppc_hc;

    // Reduced temperatures, alpha function, m parameter
    let mut trs = [0.0; 5];
    let mut alpha = [0.0; 5];
    for i in 0..5 {
        trs[i] = deg_r / tcs[i];
        let m = 0.37464 + 1.54226 * BNS_ACF[i] - 0.26992 * BNS_ACF[i] * BNS_ACF[i];
        let sqrt_tr = trs[i].sqrt();
        alpha[i] = (1.0 + m * (1.0 - sqrt_tr)).powi(2);
    }

    // BIP matrix
    let kij = calc_bips(deg_r, tpc_hc);

    // EOS parameters for this pressure
    let mut ai = [0.0; 5];
    let mut bi = [0.0; 5];
    for i in 0..5 {
        let pr_i = p_psia / pcs[i];
        ai[i] = BNS_OMEGAA[i] * alpha[i] * pr_i / (trs[i] * trs[i]);
        bi[i] = BNS_OMEGAB[i] * pr_i / trs[i];
    }

    // Mixing rules
    let mut sqrt_ai = [0.0; 5];
    for i in 0..5 {
        sqrt_ai[i] = ai[i].sqrt();
    }

    // A = sum_i sum_j zi*zj*sqrt(ai*aj)*(1-kij)
    let mut a_mix = 0.0;
    for i in 0..5 {
        for j in 0..5 {
            a_mix += zi[i] * zi[j] * sqrt_ai[i] * sqrt_ai[j] * (1.0 - kij[i][j]);
        }
    }

    // B = sum_i zi*bi
    let mut b_mix = 0.0;
    for i in 0..5 {
        b_mix += zi[i] * bi[i];
    }

    // Cubic coefficients: Z^3 + c2*Z^2 + c1*Z + c0 = 0
    let c2 = -(1.0 - b_mix);
    let c1 = a_mix - 3.0 * b_mix * b_mix - 2.0 * b_mix;
    let c0 = -(a_mix * b_mix - b_mix * b_mix - b_mix * b_mix * b_mix);

    // Solve cubic - get max (vapor) root
    let z_max = halley_cubic(c2, c1, c0);

    // Check for 3-root case (fugacity-based root selection)
    let p_d = (3.0 * c1 - c2 * c2) / 3.0;
    let q_d = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1 + 27.0 * c0) / 27.0;
    let disc = q_d * q_d / 4.0 + p_d * p_d * p_d / 27.0;

    let z_selected = if disc < -1e-15 {
        // 3 real roots - find min root via deflation
        let b_q = c2 + z_max;
        let c_q = c1 + z_max * b_q;
        let det = (b_q * b_q - 4.0 * c_q).max(0.0);
        let z_min = (-b_q - det.sqrt()) / 2.0;

        // Fugacity comparison (Gibbs criterion)
        if z_min > b_mix {
            let sqrt2: f64 = std::f64::consts::SQRT_2;
            let s2p1 = 1.0 + sqrt2;
            let s2m1 = sqrt2 - 1.0;

            let ln_phi = |zv: f64| -> f64 {
                (zv - 1.0) - (zv - b_mix).ln()
                    - a_mix / (2.0 * sqrt2 * b_mix)
                        * ((zv + s2p1 * b_mix) / (zv - s2m1 * b_mix)).ln()
            };

            if ln_phi(z_min) < ln_phi(z_max) {
                z_min
            } else {
                z_max
            }
        } else {
            z_max
        }
    } else {
        z_max
    };

    // Volume translation
    let mut vshift = 0.0;
    for i in 0..5 {
        vshift += zi[i] * BNS_VSHIFT[i] * bi[i];
    }

    z_selected - vshift
}

/// BNS full pipeline: BNS pseudocritical -> BNS PR-EOS Z-factor
#[pyfunction]
pub fn bns_zfactor_full(
    p_psia: f64,
    t_degf: f64,
    sg: f64,
    co2_frac: f64,
    h2s_frac: f64,
    n2_frac: f64,
    h2_frac: f64,
) -> PyResult<f64> {
    let deg_r = t_degf + DEGF2R;
    let (tpc, ppc, _sg_hc) = critical_properties::bns_pseudocritical_internal(
        sg, co2_frac, h2s_frac, n2_frac, h2_frac
    );

    Ok(bns_zfactor_core(p_psia, deg_r, co2_frac, h2s_frac, n2_frac, h2_frac, tpc, ppc))
}

// =========================================================================
// Public wrappers for cross-module access (used by pseudopressure)
// =========================================================================

pub fn dak_core_pub(pr: f64, tr: f64) -> f64 {
    dak_core(pr, tr)
}

pub fn hy_core_pub(pr: f64, tr: f64) -> f64 {
    hy_core(pr, tr)
}

pub fn bns_zfactor_core_pub(
    p_psia: f64,
    deg_r: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
    tpc_hc: f64,
    ppc_hc: f64,
) -> f64 {
    bns_zfactor_core(p_psia, deg_r, co2, h2s, n2, h2, tpc_hc, ppc_hc)
}
