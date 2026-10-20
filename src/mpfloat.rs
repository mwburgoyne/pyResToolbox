//! Arbitrary-precision float for the GWR inversion, on the pure-Rust
//! `astro-float` crate (replaces the `rug` GMP/MPFR bindings).
//!
//! `Mp` carries its precision in bits. Binary operations round to nearest-even
//! at the larger of the two operand precisions; scalar operations keep the
//! precision of the `Mp` operand. Transcendental functions and pi come from
//! astro-float through a per-thread constants cache. astro-float has no Euler
//! gamma, so `Mp::euler` computes it by Brent-McMillan and caches it.

use astro_float::{BigFloat, Consts, RoundingMode, Sign};
use std::cell::RefCell;
use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Sub};

const RM: RoundingMode = RoundingMode::ToEven;

thread_local! {
    static CONSTS: RefCell<Consts> = RefCell::new(Consts::new().expect("astro-float constants cache"));
    static EULER_CACHE: RefCell<Option<(usize, BigFloat)>> = const { RefCell::new(None) };
}

#[derive(Clone, Debug)]
pub struct Mp {
    v: BigFloat,
    p: usize,
}

impl Mp {
    pub fn from_f64(val: f64, prec: u32) -> Self {
        let p = prec as usize;
        Mp { v: BigFloat::from_f64(val, p), p }
    }

    pub fn from_u64(val: u64, prec: u32) -> Self {
        let p = prec as usize;
        Mp { v: BigFloat::from_word(val, p), p }
    }

    pub fn zero(prec: u32) -> Self {
        Self::from_u64(0, prec)
    }

    pub fn prec(&self) -> u32 {
        self.p as u32
    }

    /// The same value rounded (or exactly widened) to `prec` bits.
    pub fn with_prec(&self, prec: u32) -> Self {
        let p = prec as usize;
        let mut v = self.v.clone();
        v.set_precision(p, RM).expect("set_precision");
        Mp { v, p }
    }

    pub fn is_zero(&self) -> bool {
        self.v.is_zero()
    }

    pub fn is_finite(&self) -> bool {
        !self.v.is_nan() && !self.v.is_inf()
    }

    /// Binary exponent e with |self| = 0.1xxx (binary) * 2^e, None for zero/NaN/inf.
    pub fn exponent(&self) -> Option<i32> {
        if self.v.is_zero() { None } else { self.v.exponent() }
    }

    pub fn abs(&self) -> Self {
        if self.v.is_negative() { -self } else { self.clone() }
    }

    pub fn sqrt(&self) -> Self {
        Mp { v: self.v.sqrt(self.p, RM), p: self.p }
    }

    pub fn exp(&self) -> Self {
        let v = CONSTS.with(|cc| self.v.exp(self.p, RM, &mut cc.borrow_mut()));
        Mp { v, p: self.p }
    }

    pub fn ln(&self) -> Self {
        let v = CONSTS.with(|cc| self.v.ln(self.p, RM, &mut cc.borrow_mut()));
        Mp { v, p: self.p }
    }

    pub fn pi(prec: u32) -> Self {
        let p = prec as usize;
        let v = CONSTS.with(|cc| cc.borrow_mut().pi(p, RM));
        Mp { v, p }
    }

    /// Euler's gamma to `prec` bits (Brent-McMillan B1, cached per precision).
    pub fn euler(prec: u32) -> Self {
        let p = prec as usize;
        if let Some(v) = EULER_CACHE.with(|c| c.borrow().as_ref().filter(|(cp, _)| *cp >= p).map(|(_, v)| v.clone())) {
            return Mp { v, p: p.max(1) }.with_prec(prec);
        }
        let g = euler_brent_mcmillan(prec);
        EULER_CACHE.with(|c| *c.borrow_mut() = Some((p, g.v.clone())));
        g
    }

    /// Nearest f64 (round to 53 bits, then assemble exactly from the mantissa).
    pub fn to_f64(&self) -> f64 {
        if self.v.is_nan() {
            return f64::NAN;
        }
        if self.v.is_inf() {
            return if self.v.is_negative() { f64::NEG_INFINITY } else { f64::INFINITY };
        }
        if self.v.is_zero() {
            return 0.0;
        }
        let mut r = self.v.clone();
        r.set_precision(64, RM).expect("set_precision");
        r.set_precision(53, RM).expect("set_precision");
        let (words, _bits, sign, e, _) = r.as_raw_parts().expect("finite value");
        let top = *words.last().expect("mantissa word");
        // value = 0.top * 2^e = top * 2^(e - 64); top holds 53 significant bits, so this is exact
        let mag = ldexp((top >> 11) as f64, e - 53);
        if sign == Sign::Neg { -mag } else { mag }
    }
}

/// x * 2^n without the overflow of 2f64.powi(n) for |n| > 1023 (scales in steps).
fn ldexp(mut x: f64, mut n: i32) -> f64 {
    while n > 1000 {
        x *= 2f64.powi(1000);
        n -= 1000;
    }
    while n < -1000 {
        x *= 2f64.powi(-1000);
        n += 1000;
    }
    x * 2f64.powi(n)
}

/// gamma = U/V - ln(n), U = sum A_k, V = sum B_k with B_k = B_{k-1} n^2/k^2,
/// A_k = (A_{k-1} n^2/k + B_k)/k, A_0 = 0, B_0 = 1 (Brent & McMillan 1980,
/// algorithm B1). Error is about pi e^(-4n); the sums need about 2n log2(e)
/// extra bits for the growth of B_k, and run to k of about 3.59 n.
fn euler_brent_mcmillan(prec: u32) -> Mp {
    let n = (prec as f64 * std::f64::consts::LN_2 / 4.0).ceil() as u64 + 2;
    let wp = prec + (2.0 * n as f64 * std::f64::consts::LOG2_E).ceil() as u32 + 32;
    let n2 = Mp::from_u64(n * n, wp);
    let mut a = Mp::zero(wp);
    let mut b = Mp::from_u64(1, wp);
    let mut u = a.clone();
    let mut v = b.clone();
    let kmax = (3.6 * n as f64) as u64 + 10;
    for k in 1..=kmax {
        b = &(&b * &n2) / (k * k);
        a = &(&(&(&a * &n2) / k) + &b) / k;
        u = &u + &a;
        v = &v + &b;
    }
    let gamma = &(&u / &v) - &Mp::from_u64(n, wp).ln();
    gamma.with_prec(prec)
}

macro_rules! binop {
    ($tr:ident, $m:ident) => {
        impl $tr<&Mp> for &Mp {
            type Output = Mp;
            fn $m(self, rhs: &Mp) -> Mp {
                let p = self.p.max(rhs.p);
                Mp { v: self.v.$m(&rhs.v, p, RM), p }
            }
        }
        impl $tr<u64> for &Mp {
            type Output = Mp;
            fn $m(self, rhs: u64) -> Mp {
                Mp { v: self.v.$m(&BigFloat::from_word(rhs, 64), self.p, RM), p: self.p }
            }
        }
        impl $tr<f64> for &Mp {
            type Output = Mp;
            fn $m(self, rhs: f64) -> Mp {
                Mp { v: self.v.$m(&BigFloat::from_f64(rhs, 64), self.p, RM), p: self.p }
            }
        }
    };
}
binop!(Add, add);
binop!(Sub, sub);
binop!(Mul, mul);
binop!(Div, div);

impl Neg for &Mp {
    type Output = Mp;
    fn neg(self) -> Mp {
        Mp { v: BigFloat::neg(&self.v), p: self.p }
    }
}

impl PartialEq for Mp {
    fn eq(&self, other: &Mp) -> bool {
        self.v.cmp(&other.v) == Some(0)
    }
}

impl PartialOrd for Mp {
    fn partial_cmp(&self, other: &Mp) -> Option<Ordering> {
        self.v.cmp(&other.v).map(|c| c.cmp(&0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f64_round_trip() {
        for v in [1.0, 1.5, -2.75, 0.1, 1e-300, 6.02214076e23, f64::MAX / 3.0, 123456789.12345679] {
            assert_eq!(Mp::from_f64(v, 200).to_f64(), v, "{}", v);
        }
        // 1/3 at 200 bits rounds to the f64 nearest 1/3
        let third = &Mp::from_u64(1, 200) / 3u64;
        assert_eq!(third.to_f64(), 1.0 / 3.0);
    }

    #[test]
    fn euler_matches_reference() {
        // mpmath, mp.dps = 110
        let reference = "0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463";
        for digits in [20u32, 60, 100] {
            let prec = (digits as f64 * std::f64::consts::LOG2_10).ceil() as u32 + 10;
            let g = Mp::euler(prec);
            let r = Mp { v: CONSTS.with(|cc| BigFloat::parse(reference, astro_float::Radix::Dec, prec as usize + 64, RM, &mut cc.borrow_mut())), p: prec as usize + 64 };
            let rel = (&(&g - &r) / &r).abs();
            assert!(rel.is_zero() || -rel.to_f64().log2() > prec as f64 - 4.0, "{} digits", digits);
        }
    }

    #[test]
    fn pi_exp_ln() {
        let prec = 300;
        let one = Mp::from_u64(1, prec);
        assert!((&one.exp().ln() - &one).abs().to_f64() < 1e-85);
        assert_eq!(Mp::pi(prec).to_f64(), std::f64::consts::PI);
    }
}
