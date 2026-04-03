/// Alpha functions for the Peng-Robinson EOS.
///
/// Three alpha functions:
/// - alpha_water_soreide: S&W modified alpha for water/brine (Eq 9)
/// - alpha_water_mc3: Mathias-Copeman 3-parameter alpha for pure water
/// - alpha_standard_pr: Standard PR alpha for gases/hydrocarbons

/// Soreide-Whitson modified alpha for water/brine (Equation 9).
///
/// # Arguments
/// * `tr` - Reduced temperature T/Tc_water
/// * `salinity_molal` - Salt concentration in mol/kg water
pub fn alpha_water_soreide(tr: f64, salinity_molal: f64) -> f64 {
    let tr_safe = tr.max(0.01);
    let sqrt_alpha = 1.0
        + 0.4530 * (1.0 - tr_safe * (1.0 - 0.0103 * salinity_molal.powf(1.1)))
        + 0.0034 * (tr_safe.powi(-3) - 1.0);
    sqrt_alpha * sqrt_alpha
}

/// Mathias-Copeman 3-parameter alpha for pure water.
///
/// Refitted to IAPWS-95 vapor pressure data (Wagner & Pruss 2002) using
/// S&W critical properties (Tc=647.3 K, Pc=221.2 bar).
///
/// # Arguments
/// * `tr` - Reduced temperature T/Tc_water (Tc=647.3 K)
pub fn alpha_water_mc3(tr: f64) -> f64 {
    let x = 1.0 - tr.max(0.01).sqrt();
    let sqrt_alpha = 1.0 + 0.9110 * x - 0.2756 * x * x + 0.2185 * x * x * x;
    sqrt_alpha * sqrt_alpha
}

/// Standard Peng-Robinson alpha function.
///
/// # Arguments
/// * `tr` - Reduced temperature T/Tc
/// * `omega` - Acentric factor
pub fn alpha_standard_pr(tr: f64, omega: f64) -> f64 {
    let m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    let val = 1.0 + m * (1.0 - tr.max(0.0).sqrt());
    val * val
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alpha_water_soreide_freshwater() {
        // At Tr=0.5, m=0, should give a reasonable value
        let a = alpha_water_soreide(0.5, 0.0);
        assert!(a > 0.0);
        assert!(a > 1.0); // alpha > 1 below Tc
    }

    #[test]
    fn test_alpha_water_mc3() {
        // At Tr=1.0, alpha should be close to 1
        let a = alpha_water_mc3(1.0);
        assert!((a - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_alpha_standard_pr() {
        // At Tr=1.0, alpha should be 1.0 (exactly)
        let a = alpha_standard_pr(1.0, 0.3);
        assert!((a - 1.0).abs() < 1e-10);
    }
}
