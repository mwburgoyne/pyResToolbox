/// Rachford-Rice solver using the method of Nielsen & Lia (2022).
///
/// Gracefully handles catastrophic numerical roundoff errors through
/// a transformed variable approach.
///
/// Reference:
///     M. Nielsen & H. Lia, "Generalized Rachford-Rice Algorithm",
///     Fluid Phase Equilibria (2022)

/// Solve the Rachford-Rice equation.
///
/// # Arguments
/// * `zi` - Molar composition (will be normalized)
/// * `ki` - K-values for each component
/// * `tol` - Solution tolerance (default 1e-15)
/// * `max_iter` - Maximum iterations (default 100)
///
/// # Returns
/// (n_it, yi, xi, v, l) - iterations, vapor comp, liquid comp, vapor frac, liquid frac
pub fn rr_solver(
    zi: &[f64],
    ki: &[f64],
    tol: f64,
    max_iter: usize,
) -> (usize, Vec<f64>, Vec<f64>, f64, f64) {
    let nc = zi.len();
    assert_eq!(nc, ki.len());

    // Normalize feed compositions
    let z_sum: f64 = zi.iter().sum();
    let z: Vec<f64> = zi.iter().map(|&x| x / z_sum).collect();

    // Perturb K-values exactly equal to 1.0 to avoid singularity
    let k: Vec<f64> = ki
        .iter()
        .map(|&x| {
            if (x - 1.0).abs() < 1e-12 {
                1.0 + 1e-12
            } else {
                x
            }
        })
        .collect();

    // Assess if solution is nearest vapor or liquid
    let rr_at_half: f64 = z
        .iter()
        .zip(k.iter())
        .map(|(&zi, &ki)| zi * (ki - 1.0) / (1.0 + 0.5 * (ki - 1.0)))
        .sum();
    let near_vapor = rr_at_half > 0.0;

    // k_hat: reciprocal if near vapor
    let k_hat: Vec<f64> = if near_vapor {
        k.iter().map(|&x| 1.0 / x).collect()
    } else {
        k.clone()
    };

    // ci = 1 / (1 - k_hat_i)  (Eq 10)
    let ci: Vec<f64> = k_hat.iter().map(|&x| 1.0 / (1.0 - x)).collect();

    // Transformed variable bounds (Eq 11a, 11b)
    let k_hat_min = k_hat
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let k_hat_max = k_hat
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    let phi_max = (1.0 / (1.0 - k_hat_min)).min(0.5); // Eq 11a
    let phi_min = 1.0 / (1.0 - k_hat_max);             // Eq 11b

    let mut b_min = 1.0 / (phi_max - phi_min);           // Eq 15
    let mut b_max = f64::INFINITY;

    let mut b = 1.0 / (0.25 - phi_min);

    // h(b) = sum(zi * b / (1 + b*(phi_min - ci)))  (Eq 12b)
    let h = |b: f64| -> f64 {
        z.iter()
            .zip(ci.iter())
            .map(|(&zi, &ci)| zi * b / (1.0 + b * (phi_min - ci)))
            .sum()
    };

    // dh(b) = sum(zi / (1 + b*(phi_min - ci))^2)  (Eq 16b)
    let dh = |b: f64| -> f64 {
        z.iter()
            .zip(ci.iter())
            .map(|(&zi, &ci)| {
                let denom = 1.0 + b * (phi_min - ci);
                zi / (denom * denom)
            })
            .sum()
    };

    let mut n_it = 0;
    loop {
        n_it += 1;
        let h_b = h(b);
        let dh_b = dh(b);

        if h_b > 0.0 {
            b_max = b;
        } else {
            b_min = b;
        }

        if dh_b.abs() > 1e-30 {
            b -= h_b / dh_b;
        }

        if b < b_min || b > b_max {
            b = (b_min + b_max) / 2.0;
        }

        if h_b.abs() <= tol || n_it > max_iter {
            break;
        }
    }

    // Recover compositions from transformed variables (Eq 27b)
    let ui: Vec<f64> = z
        .iter()
        .zip(ci.iter())
        .map(|(&zi, &ci)| -zi * ci * b / (1.0 + b * (phi_min - ci)))
        .collect();

    // phi = (1 + b*phi_min) / b  (rearranged Eq 14b)
    let phi = (1.0 + b * phi_min) / b;

    let (v, l, yi, xi);
    if near_vapor {
        l = phi;
        v = 1.0 - l;
        yi = ui.clone();
        xi = k_hat
            .iter()
            .zip(ui.iter())
            .map(|(&kh, &u)| kh * u)
            .collect();
    } else {
        v = phi;
        l = 1.0 - v;
        xi = ui.clone();
        yi = k_hat
            .iter()
            .zip(ui.iter())
            .map(|(&kh, &u)| kh * u)
            .collect();
    }

    (n_it, yi, xi, v, l)
}

/// Convenience wrapper around rr_solver for flash calculations.
///
/// Handles single-phase detection before calling the robust solver.
///
/// # Returns
/// (V, x, y) - vapor fraction, liquid mole fractions, vapor mole fractions
pub fn solve_rachford_rice(z: &[f64], k: &[f64]) -> (f64, Vec<f64>, Vec<f64>) {
    let nc = z.len();
    assert_eq!(nc, k.len());

    // Normalize z
    let z_sum: f64 = z.iter().sum();
    let z_norm: Vec<f64> = z.iter().map(|&x| x / z_sum).collect();

    // Single-phase checks
    let sum_z_km1: f64 = z_norm
        .iter()
        .zip(k.iter())
        .map(|(&zi, &ki)| zi * (ki - 1.0))
        .sum();

    if sum_z_km1 <= 0.0 {
        // All liquid — V = 0
        let k_z: Vec<f64> = z_norm.iter().zip(k.iter()).map(|(&zi, &ki)| ki * zi).collect();
        let k_z_sum: f64 = k_z.iter().sum();
        let y: Vec<f64> = k_z.iter().map(|&x| x / k_z_sum).collect();
        return (0.0, z_norm, y);
    }

    let sum_z_km1_over_k: f64 = z_norm
        .iter()
        .zip(k.iter())
        .map(|(&zi, &ki)| zi * (ki - 1.0) / ki)
        .sum();

    if sum_z_km1_over_k >= 0.0 {
        // All vapor — V = 1
        let z_over_k: Vec<f64> = z_norm.iter().zip(k.iter()).map(|(&zi, &ki)| zi / ki).collect();
        let z_over_k_sum: f64 = z_over_k.iter().sum();
        let x: Vec<f64> = z_over_k.iter().map(|&x| x / z_over_k_sum).collect();
        return (1.0, x, z_norm);
    }

    // Two-phase: use robust solver
    let (_n_it, yi, xi, v, _l) = rr_solver(&z_norm, k, 1e-15, 100);
    (v, xi, yi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rr_solver_basic() {
        // 4-component test from Python code
        let z = [0.90, 0.05, 0.03, 0.02];
        let k = [0.001, 50.0, 80.0, 120.0];
        let (v, x, y) = solve_rachford_rice(&z, &k);

        // V should be small (mostly liquid)
        assert!(v > 0.0 && v < 1.0);
        assert_eq!(x.len(), 4);
        assert_eq!(y.len(), 4);

        // Mass balance check: z = V*y + (1-V)*x for each component
        for i in 0..4 {
            let z_sum_norm: f64 = z.iter().sum();
            let z_i = z[i] / z_sum_norm;
            let recon = v * y[i] + (1.0 - v) * x[i];
            assert!((z_i - recon).abs() < 1e-8, "Mass balance failed for component {}", i);
        }
    }

    #[test]
    fn test_all_liquid() {
        // All K < 1 should give V = 0
        let z = [0.5, 0.5];
        let k = [0.1, 0.2];
        let (v, _x, _y) = solve_rachford_rice(&z, &k);
        assert!((v - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_all_vapor() {
        // All K > 1 with sum z*(K-1)/K >= 0 should give V = 1
        let z = [0.5, 0.5];
        let k = [10.0, 20.0];
        let (v, _x, _y) = solve_rachford_rice(&z, &k);
        assert!((v - 1.0).abs() < 1e-10);
    }
}
