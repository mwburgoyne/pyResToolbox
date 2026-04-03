/// Fugacity coefficient calculations for the PR EOS.
///
/// Implements the standard PR fugacity coefficient formula and the fast
/// vectorized version used in flash_tp iterations.

use crate::vle::eos::solve_cubic_eos;

const SQRT2: f64 = 1.4142135623730951;

/// Calculate fugacity coefficient for component i in mixture.
///
/// # Arguments
/// * `z` - Compressibility factor
/// * `a_mix` - Mixture A parameter
/// * `b_mix` - Mixture B parameter
/// * `bi_over_b` - Bi/B ratio for component i
/// * `sum_xa_over_a` - (sum_j x_j * Aij) / A for component i
pub fn calc_fugacity_coeff(
    z: f64,
    a_mix: f64,
    b_mix: f64,
    bi_over_b: f64,
    sum_xa_over_a: f64,
) -> f64 {
    if b_mix < 1e-15 || a_mix < 1e-15 || z <= b_mix {
        return 1.0;
    }

    let term1 = bi_over_b * (z - 1.0);
    let term2 = -((z - b_mix).max(1e-15)).ln();

    let log_arg = (z + (1.0 + SQRT2) * b_mix) / (z + (1.0 - SQRT2) * b_mix);
    let term3 = if log_arg > 0.0 {
        (a_mix / (2.0 * SQRT2 * b_mix)) * (bi_over_b - 2.0 * sum_xa_over_a) * log_arg.ln()
    } else {
        0.0
    };

    let ln_phi = term1 + term2 + term3;
    (ln_phi.min(50.0)).exp()
}

/// Fast fugacity coefficient calculation with precomputed T/P-dependent quantities.
///
/// Used by flash_tp to avoid recomputing ai, bi, Ai, Bi on every iteration
/// when only compositions change.
///
/// # Arguments
/// * `comp` - Mole fractions (length nc)
/// * `ai_dim` - Dimensionless Ai = ai*P/(RT)^2 (length nc)
/// * `bi_dim` - Dimensionless Bi = bi*P/(RT) (length nc)
/// * `sqrt_ai` - sqrt(Ai) (length nc)
/// * `onemk` - Flattened row-major nc x nc matrix of (1 - kij)
/// * `nc` - Number of components
/// * `phase` - "liquid" or "vapor"
///
/// Returns vector of fugacity coefficients (length nc).
pub fn calc_fugacity_fast(
    comp: &[f64],
    _ai_dim: &[f64],
    bi_dim: &[f64],
    sqrt_ai: &[f64],
    onemk: &[f64],
    nc: usize,
    liquid: bool,
) -> Vec<f64> {
    // Build Aij = outer(sqrt_Ai, sqrt_Ai) * (1 - kij)
    // Then compute mixing rules
    let mut a_mix = 0.0;
    let mut b_mix = 0.0;

    // sum_xa[i] = sum_j comp[j] * Aij[i][j]
    let mut sum_xa = vec![0.0; nc];

    for i in 0..nc {
        b_mix += comp[i] * bi_dim[i];
        for j in 0..nc {
            let aij = sqrt_ai[i] * sqrt_ai[j] * onemk[i * nc + j];
            a_mix += comp[i] * comp[j] * aij;
            sum_xa[i] += comp[j] * aij;
        }
    }

    if b_mix < 1e-15 || a_mix < 1e-15 {
        return vec![1.0; nc];
    }

    let roots = solve_cubic_eos(a_mix, b_mix);
    let z = if liquid { roots[0] } else { roots[roots.len() - 1] };

    let mut phi = Vec::with_capacity(nc);
    for i in 0..nc {
        let bi_over_b = bi_dim[i] / b_mix;
        let sum_xa_over_a = sum_xa[i] / a_mix;
        phi.push(calc_fugacity_coeff(z, a_mix, b_mix, bi_over_b, sum_xa_over_a));
    }
    phi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fugacity_coeff_ideal() {
        // At very low A, B should give phi near 1
        let phi = calc_fugacity_coeff(0.999, 1e-6, 1e-7, 1.0, 0.5);
        assert!((phi - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_fugacity_fast_two_component() {
        // Simple 2-component test with no BIPs
        let comp = [0.5, 0.5];
        let ai_dim = [0.1, 0.2];
        let bi_dim = [0.01, 0.02];
        let sqrt_ai = [0.1_f64.sqrt(), 0.2_f64.sqrt()];
        let onemk = [1.0, 1.0, 1.0, 1.0]; // no BIP
        let phi = calc_fugacity_fast(&comp, &ai_dim, &bi_dim, &sqrt_ai, &onemk, 2, true);
        assert_eq!(phi.len(), 2);
        assert!(phi[0] > 0.0);
        assert!(phi[1] > 0.0);
    }
}
