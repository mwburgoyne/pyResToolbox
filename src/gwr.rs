//! Gaver-Wynn-Rho (GWR) Inverse Laplace Transform in arbitrary precision.
//!
//! Implements the GWR algorithm with arbitrary-precision arithmetic via
//! `mpfloat::Mp` (pure-Rust `astro-float`). This is the gold-standard ILT
//! method for reservoir engineering, accurate even for oscillatory transforms.
//!
//! Includes arbitrary-precision modified Bessel functions (I_0, I_1, K_0, K_1)
//! using power series for small arguments and asymptotic expansions for large,
//! with exponential scaling to avoid overflow.
//!
//! References:
//!   Valko & Abate (2004), Computers and Mathematics with Application 48(3).

use pyo3::prelude::*;
use crate::mpfloat::Mp;
use rayon::prelude::*;

// =========================================================================
//  Precision helpers
// =========================================================================

#[inline]
fn dec_to_bits(dec_digits: u32) -> u32 {
    (dec_digits as f64 * std::f64::consts::LOG2_10).ceil() as u32 + 10
}

// =========================================================================
//  Arbitrary-precision modified Bessel functions, exponentially scaled
// =========================================================================
//
// We compute scaled forms to avoid overflow for large arguments:
//   I_ne(x) = I_n(x) * exp(-x)
//   K_ne(x) = K_n(x) * exp(x)
//
// Small x: power series for I_n(x), then multiply by exp(-x) / exp(x).
// Large x: asymptotic expansion gives scaled forms directly.
//
// The switch point depends on the working precision. The asymptotic series,
// optimally truncated, is only good to a relative e^(-2x), so it is used only
// once that is below 2^-prec (with a margin). Below the switch the K_n series
// cancels terms of size e^x down to a result of size e^-x, so the series runs
// with 2x*log2(e) guard bits. A fixed x > 25 switch with x*log2(e) guard bits
// limited the transform to about 22 digits; GWR needs about 2.1*M digits.

const ASYMP_MARGIN_BITS: f64 = 16.0;
const SERIES_EXTRA_GUARD_BITS: u32 = 20;

/// True when the asymptotic expansion reaches `prec` bits at `x`.
fn use_asymptotic(x: &Mp, prec: u32) -> bool {
    2.0 * x.to_f64() * std::f64::consts::LOG2_E > prec as f64 + ASYMP_MARGIN_BITS
}

/// True when `term` no longer changes `sum` at `prec` bits.
fn series_negligible(term: &Mp, sum: &Mp, prec: u32) -> bool {
    match (term.exponent(), sum.exponent()) {
        (Some(te), Some(se)) => te < se - prec as i32,
        _ => term.is_zero(),
    }
}

/// Loop bound for an asymptotic sum: past the optimal truncation index (~2x).
/// A fixed 200 terms fell short once the precision-scaled switch passed x = 100
/// (above about 300 bits, i.e. M of 45 or more).
fn asymptotic_max_terms(x: &Mp) -> u64 {
    (4.0 * x.to_f64()) as u64 + 100
}

#[cfg(test)]
/// Scaled Bessel value from a power series: series(x) * exp(exp_sign * x),
/// evaluated with enough guard bits to absorb the K_n cancellation.
fn scaled_series(x: &Mp, prec: u32, series: fn(&Mp, u32) -> Mp, exp_sign: i32) -> Mp {
    let guard = (2.0 * x.to_f64() * std::f64::consts::LOG2_E).ceil() as u32 + SERIES_EXTRA_GUARD_BITS;
    let hp = prec + guard;
    let xhp = x.with_prec(hp);
    let val = series(&xhp, hp);
    let scale = if exp_sign < 0 { (-&xhp).exp() } else { xhp.exp() };
    (&val * &scale).with_prec(prec)
}

// --- Asymptotic expansion for large x ---
// I_nu(x) ~ e^x / sqrt(2 pi x) * sum_k (-1)^k a_k(nu) / x^k   (A&S 9.7.1)
// K_nu(x) ~ sqrt(pi/(2x)) * e^(-x) * sum_k a_k(nu) / x^k      (A&S 9.7.2)
// a_k(nu) = prod_{j=1}^k (4 nu^2 - (2j-1)^2) / (k! 8^k)
//
// Scaled forms (e^x cancels):
// I_ne(x) = 1/sqrt(2 pi x) * sum_k (-1)^k a_k / x^k
// K_ne(x) = sqrt(pi/(2x)) * sum_k a_k / x^k
//
// `asymptotic_sum` folds the factor of each term into `factor(k)`:
//   I_0: +(2k-1)^2   I_1: (2k-1)^2 - 4   (the (-1)^k of the I expansion)
//   K_0: -(2k-1)^2   K_1: 4 - (2k-1)^2

/// sum_k a_k / x^k with a_k = a_{k-1} * factor(k) / (8k), optimally truncated.
fn asymptotic_sum(x: &Mp, prec: u32, factor: fn(i64) -> i64) -> Mp {
    let one_over_x = &Mp::from_u64(1, prec) / x;
    let mut sum = Mp::from_u64(1, prec);
    let mut a_k = Mp::from_u64(1, prec);
    let mut prev_term_abs: Option<Mp> = None;

    for k in 1u64..asymptotic_max_terms(x) {
        let sq = ((2 * k - 1) * (2 * k - 1)) as i64;
        a_k = &(&(&a_k * factor(sq) as f64) / (8 * k)) * &one_over_x;

        let term_abs = a_k.abs();
        let growing = prev_term_abs.as_ref().is_some_and(|p| term_abs > *p);
        if growing || series_negligible(&a_k, &sum, prec) {
            break;
        }
        prev_term_abs = Some(term_abs);
        sum = &sum + &a_k;
    }
    sum
}

fn i_prefactor(x: &Mp, prec: u32) -> Mp {
    // 1/sqrt(2 pi x)
    &Mp::from_u64(1, prec) / &(&(x * &Mp::pi(prec)) * 2u64).sqrt()
}

fn k_prefactor(x: &Mp, prec: u32) -> Mp {
    // sqrt(pi/(2x))
    (&(&Mp::pi(prec) / x) / 2u64).sqrt()
}

fn mpfr_besseli0e_asymp(x: &Mp, prec: u32) -> Mp {
    &asymptotic_sum(x, prec, |sq| sq) * &i_prefactor(x, prec)
}

fn mpfr_besseli1e_asymp(x: &Mp, prec: u32) -> Mp {
    // 1 - 3/(8x) - 15/(128x^2) - ..., A&S 9.7.1
    &asymptotic_sum(x, prec, |sq| sq - 4) * &i_prefactor(x, prec)
}

fn mpfr_besselk0e_asymp(x: &Mp, prec: u32) -> Mp {
    &asymptotic_sum(x, prec, |sq| -sq) * &k_prefactor(x, prec)
}

fn mpfr_besselk1e_asymp(x: &Mp, prec: u32) -> Mp {
    &asymptotic_sum(x, prec, |sq| 4 - sq) * &k_prefactor(x, prec)
}

// --- Power series for small x ---

/// I_0(x) via power series: sum_{k=0}^inf (x^2/4)^k / (k!)^2
fn mpfr_besseli0_series(x: &Mp, prec: u32) -> Mp {
    let x2_over4 = &(x * x) / 4u64;
    let mut result = Mp::from_u64(1, prec);
    let mut term = Mp::from_u64(1, prec);

    for k in 1u64..100_000 {
        term = &(&term * &x2_over4) / (k * k);
        result = &result + &term;
        if k > 10 && series_negligible(&term, &result, prec) {
            break;
        }
    }
    result
}

/// I_1(x) via power series: (x/2) * sum_{k=0}^inf (x^2/4)^k / (k!*(k+1)!)
fn mpfr_besseli1_series(x: &Mp, prec: u32) -> Mp {
    let x2_over4 = &(x * x) / 4u64;
    let mut sum = Mp::from_u64(1, prec);
    let mut term = Mp::from_u64(1, prec);

    for k in 1u64..100_000 {
        term = &(&term * &x2_over4) / (k * (k + 1));
        sum = &sum + &term;
        if k > 10 && series_negligible(&term, &sum, prec) {
            break;
        }
    }
    &(x / 2u64) * &sum
}

/// K_0(x) via power series (A&S 9.6.13):
/// K_0(x) = -(ln(x/2) + gamma) * I_0(x) + sum_{k=1}^inf H_k * (x^2/4)^k / (k!)^2
#[cfg(test)]
fn mpfr_besselk0_series(x: &Mp, prec: u32) -> Mp {
    besselk0_series_from_i0(x, prec, &mpfr_besseli0_series(x, prec))
}

fn besselk0_series_from_i0(x: &Mp, prec: u32, i0: &Mp) -> Mp {
    let gamma = Mp::euler(prec);
    let ln_x_over_2 = (x / 2u64).ln();

    let mut result = -&(&(&ln_x_over_2 + &gamma) * i0);

    let x2_over4 = &(x * x) / 4u64;
    let mut term = x2_over4.clone(); // k=1
    let mut h_k = Mp::from_u64(1, prec); // H_1
    result = &result + &(&h_k * &term);

    for k in 2u64..100_000 {
        term = &(&term * &x2_over4) / (k * k);
        h_k = &h_k + &(&Mp::from_u64(1, prec) / k);
        let contribution = &h_k * &term;
        result = &result + &contribution;
        if k > 10 && series_negligible(&contribution, &result, prec) {
            break;
        }
    }
    result
}

#[cfg(test)]
/// K_1(x) via Wronskian: K_1 = (1/x - I_1*K_0) / I_0
fn mpfr_besselk1_series(x: &Mp, prec: u32) -> Mp {
    let i0 = mpfr_besseli0_series(x, prec);
    let i1 = mpfr_besseli1_series(x, prec);
    let k0 = mpfr_besselk0_series(x, prec);
    let one_over_x = &Mp::from_u64(1, prec) / x;
    &(&one_over_x - &(&i1 * &k0)) / &i0
}

// --- All four scaled functions at one argument ---

/// I_0e, I_1e, K_0e, K_1e at one argument. The VEH transform needs all four at
/// x and I_1e/K_1e at x*ReD; computing them together shares the I_0, I_1 and
/// K_0 series that the K_0 and K_1 (Wronskian) series would otherwise repeat.
struct ScaledBessel {
    i0e: Mp,
    i1e: Mp,
    k0e: Mp,
    k1e: Mp,
}

fn scaled_bessel_all(x: &Mp, prec: u32) -> ScaledBessel {
    if use_asymptotic(x, prec) {
        return ScaledBessel {
            i0e: mpfr_besseli0e_asymp(x, prec),
            i1e: mpfr_besseli1e_asymp(x, prec),
            k0e: mpfr_besselk0e_asymp(x, prec),
            k1e: mpfr_besselk1e_asymp(x, prec),
        };
    }
    let (i0, i1, k0, k1, em, ep) = unscaled_series_all(x, prec);
    ScaledBessel {
        i0e: (&i0 * &em).with_prec(prec),
        i1e: (&i1 * &em).with_prec(prec),
        k0e: (&k0 * &ep).with_prec(prec),
        k1e: (&k1 * &ep).with_prec(prec),
    }
}

/// (I_1e, K_1e) only: the x*ReD side of the transform. Same values as
/// `scaled_bessel_all`; skips the two order-0 asymptotic sums it would discard
/// (below the switch K_1 needs I_0, I_1 and K_0 anyway).
fn scaled_bessel_order1(x: &Mp, prec: u32) -> (Mp, Mp) {
    if use_asymptotic(x, prec) {
        return (mpfr_besseli1e_asymp(x, prec), mpfr_besselk1e_asymp(x, prec));
    }
    let (_, i1, _, k1, em, ep) = unscaled_series_all(x, prec);
    ((&i1 * &em).with_prec(prec), (&k1 * &ep).with_prec(prec))
}

/// I_0, I_1, K_0, K_1 by power series with K_n guard bits, plus exp(-x) and
/// exp(x), all at the guarded precision.
fn unscaled_series_all(x: &Mp, prec: u32) -> (Mp, Mp, Mp, Mp, Mp, Mp) {
    // Guard bits as in scaled_series: the K_n cancellation sets them.
    let guard = (2.0 * x.to_f64() * std::f64::consts::LOG2_E).ceil() as u32 + SERIES_EXTRA_GUARD_BITS;
    let hp = prec + guard;
    let xhp = x.with_prec(hp);
    let i0 = mpfr_besseli0_series(&xhp, hp);
    let i1 = mpfr_besseli1_series(&xhp, hp);
    let k0 = besselk0_series_from_i0(&xhp, hp, &i0);
    let k1 = &(&(&Mp::from_u64(1, hp) / &xhp) - &(&i1 * &k0)) / &i0;
    let em = (-&xhp).exp();
    let ep = xhp.exp();
    (i0, i1, k0, k1, em, ep)
}

#[cfg(test)]
fn mpfr_besseli0e(x: &Mp, prec: u32) -> Mp {
    if use_asymptotic(x, prec) { mpfr_besseli0e_asymp(x, prec) } else { scaled_series(x, prec, mpfr_besseli0_series, -1) }
}

#[cfg(test)]
fn mpfr_besselk0e(x: &Mp, prec: u32) -> Mp {
    if use_asymptotic(x, prec) { mpfr_besselk0e_asymp(x, prec) } else { scaled_series(x, prec, mpfr_besselk0_series, 1) }
}

#[cfg(test)]
fn mpfr_besselk1e(x: &Mp, prec: u32) -> Mp {
    if use_asymptotic(x, prec) { mpfr_besselk1e_asymp(x, prec) } else { scaled_series(x, prec, mpfr_besselk1_series, 1) }
}

// =========================================================================
//  Arbitrary-precision VEH Laplace-domain function (exponentially scaled)
// =========================================================================

/// Laplace-domain VEH radial flow in arbitrary precision.
///
/// Uses exponentially-scaled Bessel functions to handle large ReD.
///   I_n(x) = I_ne(x) * exp(x),  K_n(x) = K_ne(x) * exp(-x)
///   Factor out exp(delta) where delta = x*(ReD-1):
///   num = K1e(xR)*I0e(x)*exp(-2*delta) + I1e(xR)*K0e(x)
///   den = s*x * [I1e(xR)*K1e(x) - K1e(xR)*I1e(x)*exp(-2*delta)]
///
/// `x` = sqrt(s) and `bx` = its Bessel values, shared across ReD by the caller.
fn laplace_pd_mpfr(s: &Mp, x: &Mp, bx: &ScaledBessel, re_d: &Mp, prec: u32) -> Mp {
    let zero = Mp::zero(prec);
    if *s <= zero {
        return zero;
    }

    let xr = x * re_d;
    let (i1e_r, k1e_r) = scaled_bessel_order1(&xr, prec);

    let delta = x * &(re_d - 1u64);
    let two_delta = &delta * 2u64;
    let exp_neg2d = if two_delta.to_f64() > 700.0 {
        Mp::zero(prec) // effectively 0
    } else {
        (-&two_delta).exp()
    };

    let num = &(&(&k1e_r * &bx.i0e) * &exp_neg2d) + &(&i1e_r * &bx.k0e);
    let inner = &(&i1e_r * &bx.k1e) - &(&(&k1e_r * &bx.i1e) * &exp_neg2d);
    let den = &(s * x) * &inner;

    if den.is_zero() {
        return zero;
    }
    &num / &den
}

// =========================================================================
//  GWR coefficient pre-computation
// =========================================================================

struct GwrCoefficients {
    gaver_coeffs: Vec<Mp>,
    binom_table: Vec<Vec<Mp>>,
    prec_bits: u32,
}

impl GwrCoefficients {
    fn new(m: usize, prec_dec: u32) -> Self {
        let prec = dec_to_bits(prec_dec);

        let max_k = 2 * m;
        let mut fac = Vec::with_capacity(max_k + 1);
        fac.push(Mp::from_u64(1, prec));
        for k in 1..=max_k {
            let next = &fac[k - 1] * k as u64;
            fac.push(next);
        }

        let mut gaver_coeffs = Vec::with_capacity(m);
        for n in 1..=m {
            let den = &(&fac[n - 1] * &fac[n - 1]) * n as u64;
            gaver_coeffs.push(&fac[2 * n] / &den);
        }

        let mut binom_table = Vec::with_capacity(m);
        for n in 1..=m {
            let mut row = Vec::with_capacity(n + 1);
            row.push(Mp::from_u64(1, prec));
            for i in 1..=n {
                let val = &(&row[i - 1] * (n - i + 1) as u64) / i as u64;
                row.push(val);
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
    fni: &[Mp],
    m: usize,
    tau: &Mp,
    coeffs: &GwrCoefficients,
) -> f64 {
    let prec = coeffs.prec_bits;

    let m1 = m;
    let mut g0: Vec<Mp> = Vec::with_capacity(m);

    for n in 1..=m {
        let binom_row = &coeffs.binom_table[n - 1];
        let mut s = Mp::zero(prec);
        for (i, binom) in binom_row.iter().enumerate().take(n + 1) {
            let idx = n + i - 1;
            let term = binom * &fni[idx];
            s = if i & 1 == 1 { &s - &term } else { &s + &term };
        }
        g0.push(&(tau * &coeffs.gaver_coeffs[n - 1]) * &s);
    }

    let mut best = g0[m1 - 1].clone();
    let mut gm: Vec<Mp> = vec![Mp::zero(prec); m1];
    let mut gp: Vec<Mp> = vec![Mp::zero(prec); m1];

    let mut broken = false;
    for k in 0..(m1 - 1) {
        for n in (0..=(m1 - 2 - k)).rev() {
            let diff = &g0[n + 1] - &g0[n];
            if diff.is_zero() {
                broken = true;
                break;
            }
            let ratio = &Mp::from_f64(k as f64 + 1.0, prec) / &diff;
            gp[n] = &gm[n + 1] + &ratio;
            if k % 2 == 1 && n == m1 - 2 - k {
                best = gp[n].clone();
            }
        }

        if broken {
            break;
        }

        let len = m1 - k;
        gm[..len].clone_from_slice(&g0[..len]);
        g0[..len].clone_from_slice(&gp[..len]);
    }

    best.to_f64()
}

// =========================================================================
//  PyO3 exports
// =========================================================================

#[pyfunction]
#[pyo3(signature = (td_array, red_array, m=7))]
pub fn influence_tables_rust(
    py: Python<'_>,
    td_array: Vec<f64>,
    red_array: Vec<f64>,
    m: usize,
) -> PyResult<Vec<Vec<f64>>> {
    let prec_dec = (2.1 * m as f64).round() as u32;
    let prec = dec_to_bits(prec_dec);
    let coeffs = GwrCoefficients::new(m, prec_dec);
    let ln2 = Mp::from_u64(2, prec).ln();
    let red_mp: Vec<Mp> = red_array.iter().map(|&re_d| Mp::from_f64(re_d, prec)).collect();
    let taus: Vec<Mp> = td_array.iter().map(|&td| &ln2 / td).collect();
    let n_k = 2 * m;
    let n_red = red_mp.len();

    // Three stages, each parallel over all of its independent items so the
    // threads stay balanced (one task per tD column left most of them idle):
    //   1. s, sqrt(s) and the Bessel values at sqrt(s), per (tD, k); shared by every ReD
    //   2. the transform, per (tD, ReD, k)
    //   3. the GWR acceleration, per (tD, ReD)
    // Every item is computed exactly as in sequential order, so the output is
    // bit-identical whatever the thread count.
    let results = py.allow_threads(|| {
        let nodes: Vec<(Mp, Mp, ScaledBessel)> = (0..taus.len() * n_k)
            .into_par_iter()
            .map(|i| {
                let s = &taus[i / n_k] * (i % n_k + 1) as u64;
                let x = s.sqrt();
                let bx = scaled_bessel_all(&x, prec);
                (s, x, bx)
            })
            .collect();
        let fni: Vec<Mp> = (0..taus.len() * n_red * n_k)
            .into_par_iter()
            .map(|i| {
                let (t, r, k) = (i / (n_red * n_k), (i / n_k) % n_red, i % n_k);
                let (s, x, bx) = &nodes[t * n_k + k];
                laplace_pd_mpfr(s, x, bx, &red_mp[r], prec)
            })
            .collect();
        let flat: Vec<f64> = (0..taus.len() * n_red)
            .into_par_iter()
            .map(|i| gwr_single(&fni[i * n_k..(i + 1) * n_k], m, &taus[i / n_red], &coeffs))
            .collect();
        // flat is tD-major; the caller expects results[ReD][tD]
        (0..n_red)
            .map(|r| (0..taus.len()).map(|t| flat[t * n_red + r]).collect())
            .collect::<Vec<Vec<f64>>>()
    });
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
    let ln2 = Mp::from_u64(2, prec_bits).ln();

    let mut results = Vec::with_capacity(times.len());
    for &t in &times {
        let tau = &ln2 / t;

        let mut fni = Vec::with_capacity(2 * m);
        for k in 1..=(2 * m) {
            let s_val = (&tau * k as u64).to_f64();
            let fs: f64 = fn_obj.call1((s_val,))?.extract()?;
            fni.push(Mp::from_f64(fs, prec_bits));
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
        let x = Mp::from_f64(1.0, prec);
        let i0e = mpfr_besseli0e(&x, prec).to_f64();
        // I_0(1)*exp(-1) = 1.2660658777520082 * 0.367879441 = 0.46575961
        let expected = 1.2660658777520082 * (-1.0_f64).exp();
        assert!((i0e - expected).abs() < 1e-8, "I_0e(1) = {}, expected {}", i0e, expected);
    }

    #[test]
    fn test_besselk0e_small() {
        let prec = 128;
        let x = Mp::from_f64(1.0, prec);
        let k0e = mpfr_besselk0e(&x, prec).to_f64();
        // K_0(1)*exp(1) = 0.42102443824070833 * 2.71828 = 1.14446
        let expected = 0.42102443824070833 * 1.0_f64.exp();
        assert!((k0e - expected).abs() < 1e-5, "K_0e(1) = {}, expected {}", k0e, expected);
    }

    #[test]
    fn test_besselk1e_small() {
        let prec = 128;
        let x = Mp::from_f64(1.0, prec);
        let k1e = mpfr_besselk1e(&x, prec).to_f64();
        let expected = 0.6019072301972346 * 1.0_f64.exp();
        assert!((k1e - expected).abs() < 1e-5, "K_1e(1) = {}, expected {}", k1e, expected);
    }

    #[test]
    fn test_bessel_scaled_large() {
        let prec = 128;
        let x = Mp::from_f64(100.0, prec);
        let i0e = mpfr_besseli0e(&x, prec);
        assert!(i0e.is_finite());
        assert!(i0e.to_f64() > 0.0);
        let k0e = mpfr_besselk0e(&x, prec);
        assert!(k0e.is_finite());
        assert!(k0e.to_f64() > 0.0);
    }

    /// Series and asymptotic paths agree to the working precision just past
    /// the switch point, at the precisions GWR uses for M = 12 and M = 32.
    #[test]
    fn test_bessel_series_asymptotic_overlap() {
        type Pair = (fn(&Mp, u32) -> Mp, fn(&Mp, u32) -> Mp, i32);
        let pairs: [Pair; 4] = [
            (mpfr_besseli0_series, mpfr_besseli0e_asymp, -1),
            (mpfr_besseli1_series, mpfr_besseli1e_asymp, -1),
            (mpfr_besselk0_series, mpfr_besselk0e_asymp, 1),
            (mpfr_besselk1_series, mpfr_besselk1e_asymp, 1),
        ];
        for prec in [dec_to_bits(25), dec_to_bits(67), dec_to_bits(134), dec_to_bits(300)] {
            let x_switch = (prec as f64 + ASYMP_MARGIN_BITS) / (2.0 * std::f64::consts::LOG2_E);
            let x = Mp::from_f64(x_switch + 0.5, prec);
            for (series, asymp, sign) in pairs {
                let a = scaled_series(&x, prec, series, sign);
                let b = asymp(&x, prec);
                let rel = (&(&a - &b) / &b).abs();
                let bits = -rel.to_f64().log2();
                assert!(bits > prec as f64 - 12.0, "prec {}: only {:.0} bits agree", prec, bits);
            }
        }
    }

    /// K_0 at the GWR abscissa where a 60-bit evaluation loses every digit
    /// (mpmath: K_0(23.548200450309494) = 1.52394442662604e-11).
    #[test]
    fn test_besselk0_mid_argument() {
        let prec = dec_to_bits(15);
        let x = Mp::from_f64(23.548200450309494, prec);
        let k0 = mpfr_besselk0e(&x, prec).to_f64() * (-23.548200450309494_f64).exp();
        assert!((k0 / 1.52394442662604e-11 - 1.0).abs() < 1e-13, "K_0 = {}", k0);
    }

    #[test]
    fn test_gwr_coefficients_m4() {
        let coeffs = GwrCoefficients::new(4, 20);
        assert_eq!(coeffs.gaver_coeffs.len(), 4);
        assert!((coeffs.binom_table[3][2].to_f64() - 6.0).abs() < 1e-15);
    }
}
