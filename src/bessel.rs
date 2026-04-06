//! Modified Bessel functions I_0, I_1, K_0, K_1 for f64.
//!
//! Provides both standard and exponentially-scaled versions:
//!   besseli0e(x) = I_0(x) * e^(-|x|)
//!   besselk0e(x) = K_0(x) * e^(x)
//! The scaled forms avoid overflow/underflow for large arguments.
//!
//! Polynomial approximations from Abramowitz & Stegun (1972)
//! with Cephes library coefficients.

/// I_0(x) * e^(-|x|) — exponentially scaled modified Bessel I_0
pub fn besseli0e(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let t = (x / 3.75).powi(2);
        let unscaled = 1.0 + t * (3.5156229
            + t * (3.0899424
                + t * (1.2067492
                    + t * (0.2659732
                        + t * (0.0360768
                            + t * 0.0045813)))));
        unscaled * (-ax).exp()
    } else {
        let t = 3.75 / ax;
        let ans = 0.39894228
            + t * (0.01328592
                + t * (0.00225319
                    + t * (-0.00157565
                        + t * (0.00916281
                            + t * (-0.02057706
                                + t * (0.02635537
                                    + t * (-0.01647633
                                        + t * 0.00392377)))))));
        ans / ax.sqrt()
    }
}

/// I_1(x) * e^(-|x|) — exponentially scaled modified Bessel I_1
pub fn besseli1e(x: f64) -> f64 {
    let ax = x.abs();
    let result = if ax < 3.75 {
        let t = (x / 3.75).powi(2);
        let unscaled = ax * (0.5
            + t * (0.87890594
                + t * (0.51498869
                    + t * (0.15084934
                        + t * (0.02658733
                            + t * (0.00301532
                                + t * 0.00032411))))));
        unscaled * (-ax).exp()
    } else {
        let t = 3.75 / ax;
        let ans = 0.39894228
            + t * (-0.03988024
                + t * (-0.00362018
                    + t * (0.00163801
                        + t * (-0.01031555
                            + t * (0.02282967
                                + t * (-0.02895312
                                    + t * (0.01787654
                                        + t * (-0.00420059))))))));
        ans / ax.sqrt()
    };
    if x < 0.0 { -result } else { result }
}

/// K_0(x) * e^(x) — exponentially scaled modified Bessel K_0
///
/// x must be > 0.
pub fn besselk0e(x: f64) -> f64 {
    debug_assert!(x > 0.0);
    if x <= 2.0 {
        let t = x * x / 4.0;
        let i0 = besseli0(x);
        let poly = -0.57721566
            + t * (0.42278420
                + t * (0.23069756
                    + t * (0.03488590
                        + t * (0.00262698
                            + t * (0.00010750
                                + t * 0.00000740)))));
        ((-x.ln() + std::f64::consts::LN_2 - 0.5772156649015329) * i0 + poly) * x.exp()
    } else {
        let t = 2.0 / x;
        let ans = 1.25331414
            + t * (-0.07832358
                + t * (0.02189568
                    + t * (-0.01062446
                        + t * (0.00587872
                            + t * (-0.00251540
                                + t * 0.00053208)))));
        ans / x.sqrt()
    }
}

/// K_1(x) * e^(x) — exponentially scaled modified Bessel K_1
///
/// x must be > 0.
pub fn besselk1e(x: f64) -> f64 {
    debug_assert!(x > 0.0);
    if x <= 2.0 {
        let t = x * x / 4.0;
        let i1 = besseli1(x);
        let poly = (1.0 / x)
            * (1.0
                + t * (0.15443144
                    + t * (-0.67278579
                        + t * (-0.18156897
                            + t * (-0.01919402
                                + t * (-0.00110404
                                    + t * (-0.00004686)))))));
        ((x.ln() - std::f64::consts::LN_2 + 0.5772156649015329) * i1 + poly) * x.exp()
    } else {
        let t = 2.0 / x;
        let ans = 1.25331414
            + t * (0.23498619
                + t * (-0.03655620
                    + t * (0.01504268
                        + t * (-0.00780353
                            + t * (0.00325614
                                + t * (-0.00068245))))));
        ans / x.sqrt()
    }
}

// === Unscaled versions (for small arguments) ===

/// Modified Bessel function I_0(x)
pub fn besseli0(x: f64) -> f64 {
    besseli0e(x) * x.abs().exp()
}

/// Modified Bessel function I_1(x)
pub fn besseli1(x: f64) -> f64 {
    besseli1e(x) * x.abs().exp()
}

/// Modified Bessel function K_0(x).  x must be > 0.
pub fn besselk0(x: f64) -> f64 {
    besselk0e(x) * (-x).exp()
}

/// Modified Bessel function K_1(x).  x must be > 0.
pub fn besselk1(x: f64) -> f64 {
    besselk1e(x) * (-x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_besseli0() {
        assert!((besseli0(0.0) - 1.0).abs() < 1e-10);
        assert!((besseli0(1.0) - 1.2660658777520082).abs() < 1e-7);
        assert!((besseli0(5.0) - 27.239871823604442).abs() < 1e-4);
    }

    #[test]
    fn test_besseli1() {
        assert!((besseli1(0.0)).abs() < 1e-10);
        assert!((besseli1(1.0) - 0.5651591039924851).abs() < 1e-7);
    }

    #[test]
    fn test_besselk0() {
        assert!((besselk0(1.0) - 0.42102443824070833).abs() < 1e-7);
        assert!((besselk0(5.0) - 0.0036910982234656676).abs() < 1e-10);
    }

    #[test]
    fn test_besselk1() {
        assert!((besselk1(1.0) - 0.6019072301972346).abs() < 1e-7);
        assert!((besselk1(5.0) - 0.004044613445452164).abs() < 1e-10);
    }

    #[test]
    fn test_scaled_large_arg() {
        // For x=100, I_0(100) would overflow f64, but scaled version works
        let ie0 = besseli0e(100.0);
        assert!(ie0.is_finite());
        assert!(ie0 > 0.0);

        let ke0 = besselk0e(100.0);
        assert!(ke0.is_finite());
        assert!(ke0 > 0.0);
    }
}
