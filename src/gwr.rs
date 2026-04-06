//! Gaver-Wynn-Rho (GWR) Inverse Laplace Transform — Rust + MPFR acceleration.
//!
//! Implements the GWR algorithm with arbitrary-precision arithmetic via the
//! `rug` crate (GMP/MPFR bindings).  This is the gold-standard ILT method
//! for reservoir engineering — accurate even for oscillatory transforms.
//!
//! Includes MPFR-precision modified Bessel functions (I_0, I_1, K_0, K_1)
//! using power series for small arguments and asymptotic expansions for large,
//! with exponential scaling to avoid overflow.
//!
//! References:
//!   Valko & Abate (2004), Computers and Mathematics with Application 48(3).

use pyo3::prelude::*;
use rug::Float;
use rug::float::Round;
use rug::ops::AssignRound;

// =========================================================================
//  Precision helpers
// =========================================================================

#[inline]
fn dec_to_bits(dec_digits: u32) -> u32 {
    (dec_digits as f64 * std::f64::consts::LOG2_10).ceil() as u32 + 10
}

#[inline]
fn new_float(prec: u32) -> Float {
    Float::new(prec)
}

#[inline]
fn float_from_f64(val: f64, prec: u32) -> Float {
    Float::with_val(prec, val)
}

// =========================================================================
//  MPFR Modified Bessel Functions — exponentially scaled
// =========================================================================
//
// We compute scaled forms to avoid overflow for large arguments:
//   I_ne(x) = I_n(x) * exp(-x)
//   K_ne(x) = K_n(x) * exp(x)
//
// Small x: power series for I_n(x), then multiply by exp(-x) / exp(x).
// Large x: asymptotic expansion gives scaled forms directly.
//
// Threshold: x > 20 uses asymptotic. Series converges for all x but
// K_n series has catastrophic cancellation for large x.

const BESSEL_THRESHOLD: f64 = 25.0;

// --- Asymptotic expansion for large x ---
// I_ν(x) ~ e^x / sqrt(2πx) * sum_{k=0}^K (-1)^k * a_k(ν) / x^k
// K_ν(x) ~ sqrt(π/(2x)) * e^(-x) * sum_{k=0}^K a_k(ν) / x^k
//
// a_k(ν) = prod_{j=1}^k (4ν²-(2j-1)²) / (k! * 8^k)
//
// Scaled forms (e^x cancels):
// I_ne(x) = 1/sqrt(2πx) * sum_{k=0}^K (-1)^k * a_k / x^k
// K_ne(x) = sqrt(π/(2x)) * sum_{k=0}^K a_k / x^k

fn mpfr_besseli0e_asymp(x: &Float, prec: u32) -> Float {
    // mu = 4*0^2 = 0
    let one_over_x = Float::with_val(prec, 1) / x;
    let one_over_8x = Float::with_val(prec, &one_over_x / 8u32);

    // Sum the asymptotic series: sum_{k=0}^K (-1)^k * a_k / x^k
    // a_0 = 1
    // a_k = a_{k-1} * (0-(2k-1)^2) / (k*8) = a_{k-1} * (-(2k-1)^2) / (8k)
    // For I (alternating): term_k = (-1)^k * a_k / x^k
    let mut sum = Float::with_val(prec, 1);
    let mut a_k = Float::with_val(prec, 1); // a_0 = 1
    let mut prev_term_abs = Float::with_val(prec, f64::MAX);

    for k in 1u64..200 {
        let sq = (2 * k - 1) * (2 * k - 1);
        // a_k = a_{k-1} * (-(2k-1)^2) / (8k*x)
        // But we track a_k / x^k directly:
        // factor = -(2k-1)^2 / (8*k*x)
        a_k *= Float::with_val(prec, sq);
        a_k /= 8u32 * k as u32;
        a_k *= &one_over_x;

        let term_abs = Float::with_val(prec, a_k.clone().abs());
        // Optimal truncation: stop when terms start growing
        if term_abs > prev_term_abs {
            break;
        }
        prev_term_abs.assign_round(&term_abs, Round::Nearest);

        // For I_0: alternating signs ((-1)^k) but a_k already has sign from -(2k-1)^2
        // Actually: a_k accumulates the sign from the product of -(2k-1)^2 factors
        // So a_k = prod_{j=1}^k (-(2j-1)^2) / (j*8) = (-1)^k * prod (2j-1)^2 / (j*8)
        // For I_0e, the expansion is 1/sqrt(2*pi*x) * sum (-1)^k * a_k(ν)/x^k
        // Since our a_k already has (-1)^k built in, just add it
        sum += &a_k;
    }

    // I_0e(x) = sum / sqrt(2*pi*x)
    let two_pi_x = Float::with_val(prec, x * std::f64::consts::TAU);
    Float::with_val(prec, &sum / two_pi_x.sqrt())
}

fn mpfr_besseli1e_asymp(x: &Float, prec: u32) -> Float {
    // mu = 4*1^2 = 4
    let one_over_x = Float::with_val(prec, 1) / x;

    let mut sum = Float::with_val(prec, 1);
    let mut a_k = Float::with_val(prec, 1);
    let mut prev_term_abs = Float::with_val(prec, f64::MAX);

    for k in 1u64..200 {
        let sq = (2 * k - 1) * (2 * k - 1);
        // factor for nu=1: (4 - (2k-1)^2) / (8*k*x)
        let numer = 4i64 - sq as i64;
        a_k *= Float::with_val(prec, numer);
        a_k /= 8u32 * k as u32;
        a_k *= &one_over_x;

        let term_abs = Float::with_val(prec, a_k.clone().abs());
        if term_abs > prev_term_abs {
            break;
        }
        prev_term_abs.assign_round(&term_abs, Round::Nearest);
        sum += &a_k;
    }

    let two_pi_x = Float::with_val(prec, x * std::f64::consts::TAU);
    Float::with_val(prec, &sum / two_pi_x.sqrt())
}

fn mpfr_besselk0e_asymp(x: &Float, prec: u32) -> Float {
    // mu = 0
    let one_over_x = Float::with_val(prec, 1) / x;

    let mut sum = Float::with_val(prec, 1);
    let mut a_k = Float::with_val(prec, 1);
    let mut prev_term_abs = Float::with_val(prec, f64::MAX);

    for k in 1u64..200 {
        let sq = (2 * k - 1) * (2 * k - 1);
        // For K_0: a_k(0) = prod (0-(2j-1)^2) / (j*8*x)
        // = prod (-(2j-1)^2) / (8j*x)
        // For K: no alternating sign — K_ne = sqrt(pi/(2x)) * sum a_k/x^k
        // But a_k has sign from the -(2j-1)^2 products.
        // Actually for K the formula is sum with ALL positive signs:
        // K_ne(x) = sqrt(pi/(2x)) * sum_{k=0} a_k(nu) / x^k
        // where a_k = prod_{j=1}^k (mu-(2j-1)^2) / (j! * 8^j)
        // Wait, no. The a_k are: a_0=1, a_k = prod_{j=1}^k (mu-(2j-1)^2)/(8j)
        // For mu=0: a_k = prod (-1*(2j-1)^2) / (8j) — this alternates sign
        // But K has all positive signs: sum a_k / x^k where a_k are as defined

        // Let me reconsider. From A&S 9.7.2:
        // K_0(x) ~ sqrt(pi/(2x)) * e^(-x) * (1 + sum_{k=1} mu_k / (8x)^k)
        // where mu_k = (1^2)(3^2)...(2k-1)^2 / k!  (all positive for nu=0)
        // Wait no, that's not right either.

        // Let me be precise. A&S 9.7.2:
        // K_nu(x) ~ sqrt(pi/(2x)) * e^(-x) * sum_{k=0}^inf (nu,k) / (2x)^k
        // where (nu,k) = prod_{s=0}^{k-1} (4nu^2-(2s+1)^2) / (8*(s+1))  ← this is a_k
        //
        // For nu=0: (0,k) = prod_{s=0}^{k-1} (-(2s+1)^2)/(8(s+1))
        //   k=1: -(1)/(8) = -1/8
        //   k=2: (-1/8) * (-9/16) = 9/128
        //
        // So K series has ALTERNATING signs (negative, positive, negative...)
        // But the formula says these are added: K_ne = sqrt(pi/2x) * (1 + a1 + a2 + ...)
        // where a_k alternate in sign.

        // For the implementation: a_k accumulates sign from the product
        a_k *= Float::with_val(prec, -(sq as f64));
        a_k /= 8u32 * k as u32;
        a_k *= &one_over_x;

        let term_abs = Float::with_val(prec, a_k.clone().abs());
        if term_abs > prev_term_abs {
            break;
        }
        prev_term_abs.assign_round(&term_abs, Round::Nearest);
        sum += &a_k;
    }

    let pi_over_2x = Float::with_val(prec, std::f64::consts::FRAC_PI_2) / x;
    Float::with_val(prec, &sum * pi_over_2x.sqrt())
}

fn mpfr_besselk1e_asymp(x: &Float, prec: u32) -> Float {
    // mu = 4
    let one_over_x = Float::with_val(prec, 1) / x;

    let mut sum = Float::with_val(prec, 1);
    let mut a_k = Float::with_val(prec, 1);
    let mut prev_term_abs = Float::with_val(prec, f64::MAX);

    for k in 1u64..200 {
        let sq = (2 * k - 1) * (2 * k - 1);
        let numer = 4.0 - sq as f64;
        a_k *= Float::with_val(prec, numer);
        a_k /= 8u32 * k as u32;
        a_k *= &one_over_x;

        let term_abs = Float::with_val(prec, a_k.clone().abs());
        if term_abs > prev_term_abs {
            break;
        }
        prev_term_abs.assign_round(&term_abs, Round::Nearest);
        sum += &a_k;
    }

    let pi_over_2x = Float::with_val(prec, std::f64::consts::FRAC_PI_2) / x;
    Float::with_val(prec, &sum * pi_over_2x.sqrt())
}

// --- Power series for small x ---

/// I_0(x) via power series: sum_{k=0}^inf (x^2/4)^k / (k!)^2
fn mpfr_besseli0_series(x: &Float, prec: u32) -> Float {
    let x2_over4 = Float::with_val(prec, x * x) / 4u32;
    let mut result = Float::with_val(prec, 1);
    let mut term = Float::with_val(prec, 1);

    for k in 1u64..100_000 {
        term *= &x2_over4;
        term /= k * k;
        result += &term;
        if term.is_zero() || (k > 10 && Float::with_val(prec, &term / &result).to_f64().abs() < 1e-50) {
            break;
        }
    }
    result
}

/// I_1(x) via power series: (x/2) * sum_{k=0}^inf (x^2/4)^k / (k!*(k+1)!)
fn mpfr_besseli1_series(x: &Float, prec: u32) -> Float {
    let x2_over4 = Float::with_val(prec, x * x) / 4u32;
    let mut sum = Float::with_val(prec, 1);
    let mut term = Float::with_val(prec, 1);

    for k in 1u64..100_000 {
        term *= &x2_over4;
        term /= k * (k + 1);
        sum += &term;
        if term.is_zero() || (k > 10 && Float::with_val(prec, &term / &sum).to_f64().abs() < 1e-50) {
            break;
        }
    }
    Float::with_val(prec, x / 2u32) * sum
}

/// K_0(x) via power series (A&S 9.6.13):
/// K_0(x) = -(ln(x/2) + gamma) * I_0(x) + sum_{k=1}^inf H_k * (x^2/4)^k / (k!)^2
fn mpfr_besselk0_series(x: &Float, prec: u32) -> Float {
    let i0 = mpfr_besseli0_series(x, prec);
    let gamma = Float::with_val(prec,
        Float::parse("0.57721566490153286060651209008240243104215933593992").unwrap());
    let ln_x_over_2 = Float::with_val(prec, x / 2u32).ln();

    let mut result = Float::with_val(prec, &ln_x_over_2 + &gamma);
    result = -result * &i0;

    let x2_over4 = Float::with_val(prec, x * x) / 4u32;
    let mut term = Float::with_val(prec, &x2_over4); // k=1
    let mut h_k = Float::with_val(prec, 1); // H_1
    result += Float::with_val(prec, &h_k * &term);

    for k in 2u64..100_000 {
        term *= &x2_over4;
        term /= k * k;
        h_k += Float::with_val(prec, 1) / Float::with_val(prec, k);
        let contribution = Float::with_val(prec, &h_k * &term);
        result += &contribution;
        if contribution.is_zero() || (k > 10 && Float::with_val(prec, &contribution / &result).to_f64().abs() < 1e-50) {
            break;
        }
    }
    result
}

/// K_1(x) via Wronskian: K_1 = (1/x - I_1*K_0) / I_0
fn mpfr_besselk1_series(x: &Float, prec: u32) -> Float {
    let i0 = mpfr_besseli0_series(x, prec);
    let i1 = mpfr_besseli1_series(x, prec);
    let k0 = mpfr_besselk0_series(x, prec);
    let one_over_x = Float::with_val(prec, 1) / x;
    Float::with_val(prec, &one_over_x - Float::with_val(prec, &i1 * &k0)) / &i0
}

// --- Dispatch: series or asymptotic ---

fn mpfr_besseli0e(x: &Float, prec: u32) -> Float {
    if x.to_f64() < BESSEL_THRESHOLD {
        // I_0(x)*exp(-x): I_0 ~ exp(x)/sqrt(2πx), so product ~ 1/sqrt(2πx).
        // Need guard bits ≈ x*log2(e) to preserve precision in I_0 before scaling.
        let guard = (x.to_f64() * std::f64::consts::LOG2_E).ceil() as u32 + 20;
        let hp = prec + guard;
        let xhp = Float::with_val(hp, x);
        let i0 = mpfr_besseli0_series(&xhp, hp);
        let result = Float::with_val(hp, &i0 * Float::with_val(hp, -&xhp).exp());
        Float::with_val(prec, &result)
    } else {
        mpfr_besseli0e_asymp(x, prec)
    }
}

fn mpfr_besseli1e(x: &Float, prec: u32) -> Float {
    if x.to_f64() < BESSEL_THRESHOLD {
        let guard = (x.to_f64() * std::f64::consts::LOG2_E).ceil() as u32 + 20;
        let hp = prec + guard;
        let xhp = Float::with_val(hp, x);
        let i1 = mpfr_besseli1_series(&xhp, hp);
        let result = Float::with_val(hp, &i1 * Float::with_val(hp, -&xhp).exp());
        Float::with_val(prec, &result)
    } else {
        mpfr_besseli1e_asymp(x, prec)
    }
}

fn mpfr_besselk0e(x: &Float, prec: u32) -> Float {
    if x.to_f64() < BESSEL_THRESHOLD {
        // Guard bits to compensate for cancellation in K_0 series.
        // Cancellation is O(e^x), needing ~x*log2(e) extra bits.
        let guard = (x.to_f64() * std::f64::consts::LOG2_E).ceil() as u32 + 20;
        let hp = prec + guard;
        let xhp = Float::with_val(hp, x);
        let k0 = mpfr_besselk0_series(&xhp, hp);
        let result = Float::with_val(hp, &k0 * xhp.exp());
        Float::with_val(prec, &result)
    } else {
        mpfr_besselk0e_asymp(x, prec)
    }
}

fn mpfr_besselk1e(x: &Float, prec: u32) -> Float {
    if x.to_f64() < BESSEL_THRESHOLD {
        let guard = (x.to_f64() * std::f64::consts::LOG2_E).ceil() as u32 + 20;
        let hp = prec + guard;
        let xhp = Float::with_val(hp, x);
        let k1 = mpfr_besselk1_series(&xhp, hp);
        let result = Float::with_val(hp, &k1 * xhp.exp());
        Float::with_val(prec, &result)
    } else {
        mpfr_besselk1e_asymp(x, prec)
    }
}

// =========================================================================
//  MPFR-precision VEH Laplace-domain function (exponentially scaled)
// =========================================================================

/// Laplace-domain VEH radial flow in arbitrary precision.
///
/// Uses exponentially-scaled Bessel functions to handle large ReD.
/// Same scaling approach as the f64 version:
///   I_n(x) = I_ne(x) * exp(x),  K_n(x) = K_ne(x) * exp(-x)
///   Factor out exp(delta) where delta = x*(ReD-1):
///   num = K1e(xR)*I0e(x)*exp(-2*delta) + I1e(xR)*K0e(x)
///   den = s*x * [I1e(xR)*K1e(x) - K1e(xR)*I1e(x)*exp(-2*delta)]
fn laplace_pd_mpfr(s: &Float, re_d: &Float, prec: u32) -> Float {
    let zero = new_float(prec);
    if *s <= zero {
        return zero;
    }

    let x = Float::with_val(prec, s.clone().sqrt());
    let xr = Float::with_val(prec, &x * re_d);

    let ie0x = mpfr_besseli0e(&x, prec);
    let ie1x = mpfr_besseli1e(&x, prec);
    let ke0x = mpfr_besselk0e(&x, prec);
    let ke1x = mpfr_besselk1e(&x, prec);
    let ie1xr = mpfr_besseli1e(&xr, prec);
    let ke1xr = mpfr_besselk1e(&xr, prec);

    let delta = Float::with_val(prec, &x * Float::with_val(prec, re_d - 1u32));
    let two_delta = Float::with_val(prec, &delta * 2u32);
    let exp_neg2d = if two_delta.to_f64() > 700.0 {
        new_float(prec) // effectively 0
    } else {
        Float::with_val(prec, -&two_delta).exp()
    };

    let num = Float::with_val(prec, &ke1xr * &ie0x) * &exp_neg2d
        + Float::with_val(prec, &ie1xr * &ke0x);

    let inner = Float::with_val(prec, &ie1xr * &ke1x)
        - Float::with_val(prec, &ke1xr * &ie1x) * &exp_neg2d;
    let den = Float::with_val(prec, s * &x) * inner;

    if den.is_zero() {
        return zero;
    }
    Float::with_val(prec, &num / &den)
}

// =========================================================================
//  GWR coefficient pre-computation
// =========================================================================

struct GwrCoefficients {
    gaver_coeffs: Vec<Float>,
    binom_table: Vec<Vec<Float>>,
    prec_bits: u32,
}

impl GwrCoefficients {
    fn new(m: usize, prec_dec: u32) -> Self {
        let prec = dec_to_bits(prec_dec);

        let max_k = 2 * m;
        let mut fac = Vec::with_capacity(max_k + 1);
        fac.push(Float::with_val(prec, 1));
        for k in 1..=max_k {
            let prev = &fac[k - 1];
            fac.push(Float::with_val(prec, prev * k as u32));
        }

        let mut gaver_coeffs = Vec::with_capacity(m);
        for n in 1..=m {
            let num = &fac[2 * n];
            let n_fac_m1 = &fac[n - 1];
            let den = Float::with_val(prec, n_fac_m1 * n_fac_m1) * n as u32;
            gaver_coeffs.push(Float::with_val(prec, num / &den));
        }

        let mut binom_table = Vec::with_capacity(m);
        for n in 1..=m {
            let mut row = Vec::with_capacity(n + 1);
            row.push(Float::with_val(prec, 1));
            for i in 1..=n {
                let prev = &row[i - 1];
                let val = Float::with_val(prec, prev * (n - i + 1) as u32) / i as u32;
                row.push(Float::with_val(prec, val));
            }
            binom_table.push(row);
        }

        GwrCoefficients {
            gaver_coeffs,
            binom_table,
            prec_bits: prec,
        }
    }
}

// =========================================================================
//  Core GWR algorithm
// =========================================================================

fn gwr_single(
    fni: &[Float],
    m: usize,
    tau: &Float,
    coeffs: &GwrCoefficients,
) -> f64 {
    let prec = coeffs.prec_bits;

    let m1 = m;
    let mut g0: Vec<Float> = Vec::with_capacity(m);

    for n in 1..=m {
        let binom_row = &coeffs.binom_table[n - 1];
        let mut s = new_float(prec);
        for i in 0..=n {
            let idx = n + i - 1;
            let term = Float::with_val(prec, &binom_row[i] * &fni[idx]);
            if i & 1 == 1 {
                s -= term;
            } else {
                s += term;
            }
        }
        g0.push(Float::with_val(prec, tau * &coeffs.gaver_coeffs[n - 1]) * s);
    }

    let mut best = Float::with_val(prec, &g0[m1 - 1]);
    let mut gm: Vec<Float> = vec![new_float(prec); m1];
    let mut gp: Vec<Float> = vec![new_float(prec); m1];

    let mut broken = false;
    for k in 0..(m1 - 1) {
        for n in (0..=(m1 - 2 - k)).rev() {
            let diff = Float::with_val(prec, &g0[n + 1] - &g0[n]);
            if diff.is_zero() {
                broken = true;
                break;
            }
            let ratio = Float::with_val(prec, (k as f64 + 1.0) / &diff);
            gp[n] = Float::with_val(prec, &gm[n + 1] + &ratio);
            if k % 2 == 1 && n == m1 - 2 - k {
                best.assign_round(&gp[n], Round::Nearest);
            }
        }

        if broken {
            break;
        }

        for n in 0..(m1 - k) {
            gm[n].assign_round(&g0[n], Round::Nearest);
            g0[n].assign_round(&gp[n], Round::Nearest);
        }
    }

    best.to_f64()
}

// =========================================================================
//  PyO3 exports
// =========================================================================

#[pyfunction]
#[pyo3(signature = (td_array, red_array, m=8))]
pub fn influence_tables_rust(
    td_array: Vec<f64>,
    red_array: Vec<f64>,
    m: usize,
) -> PyResult<Vec<Vec<f64>>> {
    let prec_dec = (2.1 * m as f64).round() as u32;
    let prec = dec_to_bits(prec_dec);
    let coeffs = GwrCoefficients::new(m, prec_dec);

    let results: Vec<Vec<f64>> = red_array
        .iter()
        .map(|&re_d| {
            let re_d_mpfr = float_from_f64(re_d, prec);
            td_array
                .iter()
                .map(|&td| {
                    let mut tau = Float::with_val(prec, 2u32).ln();
                    tau /= td;

                    let mut fni = Vec::with_capacity(2 * m);
                    for k in 1..=(2 * m) {
                        let s = Float::with_val(prec, &tau * k as u32);
                        fni.push(laplace_pd_mpfr(&s, &re_d_mpfr, prec));
                    }

                    gwr_single(&fni, m, &tau, &coeffs)
                })
                .collect()
        })
        .collect();

    Ok(results)
}

#[pyfunction]
#[pyo3(signature = (fn_obj, times, m=32, prec=None))]
pub fn gwr_rust(
    _py: Python<'_>,
    fn_obj: Bound<'_, PyAny>,
    times: Vec<f64>,
    m: usize,
    prec: Option<u32>,
) -> PyResult<Vec<f64>> {
    let prec_dec = prec.unwrap_or_else(|| (2.1 * m as f64).round() as u32);
    let prec_bits = dec_to_bits(prec_dec);
    let coeffs = GwrCoefficients::new(m, prec_dec);

    let mut results = Vec::with_capacity(times.len());
    for &t in &times {
        let mut tau = Float::with_val(prec_bits, 2u32).ln();
        tau /= t;

        let mut fni = Vec::with_capacity(2 * m);
        for k in 1..=(2 * m) {
            let s_val = Float::with_val(prec_bits, &tau * k as u32).to_f64();
            let fs: f64 = fn_obj.call1((s_val,))?.extract()?;
            fni.push(float_from_f64(fs, prec_bits));
        }

        results.push(gwr_single(&fni, m, &tau, &coeffs));
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_besseli0e_small() {
        let prec = 128;
        let x = Float::with_val(prec, 1.0);
        let i0e = mpfr_besseli0e(&x, prec).to_f64();
        // I_0(1)*exp(-1) = 1.2660658777520082 * 0.367879441 = 0.46575961
        let expected = 1.2660658777520082 * (-1.0_f64).exp();
        assert!((i0e - expected).abs() < 1e-8, "I_0e(1) = {}, expected {}", i0e, expected);
    }

    #[test]
    fn test_besselk0e_small() {
        let prec = 128;
        let x = Float::with_val(prec, 1.0);
        let k0e = mpfr_besselk0e(&x, prec).to_f64();
        // K_0(1)*exp(1) = 0.42102443824070833 * 2.71828 = 1.14446
        let expected = 0.42102443824070833 * 1.0_f64.exp();
        assert!((k0e - expected).abs() < 1e-5, "K_0e(1) = {}, expected {}", k0e, expected);
    }

    #[test]
    fn test_besselk1e_small() {
        let prec = 128;
        let x = Float::with_val(prec, 1.0);
        let k1e = mpfr_besselk1e(&x, prec).to_f64();
        let expected = 0.6019072301972346 * 1.0_f64.exp();
        assert!((k1e - expected).abs() < 1e-5, "K_1e(1) = {}, expected {}", k1e, expected);
    }

    #[test]
    fn test_bessel_scaled_large() {
        let prec = 128;
        let x = Float::with_val(prec, 100.0);
        let i0e = mpfr_besseli0e(&x, prec);
        assert!(i0e.is_finite());
        assert!(i0e.to_f64() > 0.0);
        let k0e = mpfr_besselk0e(&x, prec);
        assert!(k0e.is_finite());
        assert!(k0e.to_f64() > 0.0);
    }

    #[test]
    fn test_gwr_coefficients_m4() {
        let coeffs = GwrCoefficients::new(4, 20);
        assert_eq!(coeffs.gaver_coeffs.len(), 4);
        assert!((coeffs.binom_table[3][2].to_f64() - 6.0).abs() < 1e-15);
    }
}
