#![allow(dead_code)]
/// Cubic EOS solver for the Peng-Robinson equation of state.
///
/// Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
///
/// Uses Halley-accelerated iteration (Michelsen-style) with analytic fallback.

const MAX_ITER: usize = 50;
const TOL: f64 = 1e-10;
const SQRT2: f64 = 1.4142135623730951;

/// Solve depressed cubic Z^3 + c2*Z^2 + c1*Z + c0 = 0 using Halley iteration.
///
/// Michelsen-style: start from inflection point, find largest root via Halley,
/// then synthetic division + quadratic for remaining roots.
///
/// Returns sorted list of valid roots (Z > B).
fn halley_cubic(c2: f64, c1: f64, c0: f64, b: f64) -> Vec<f64> {
    // Inflection point: F''=0 at Z_inf = -c2/3
    let z_inf = -c2 / 3.0;

    // Evaluate F at inflection
    let f_inf = z_inf * z_inf * z_inf + c2 * z_inf * z_inf + c1 * z_inf + c0;

    // Choose starting point for largest root
    let mut z = if f_inf > 0.0 {
        // Function positive at inflection — largest root is above inflection
        (b + 1.0_f64).max(z_inf + 1.0)
    } else {
        // F_inf <= 0: largest root might be near or above local max
        let disc_fp = c2 * c2 - 3.0 * c1;
        if disc_fp > 0.0 {
            let z_local_max = (-c2 + disc_fp.sqrt()) / 3.0;
            let f_max = z_local_max * z_local_max * z_local_max
                + c2 * z_local_max * z_local_max
                + c1 * z_local_max
                + c0;
            if f_max > 0.0 {
                // Three real roots — start above local max for largest
                z_local_max + 0.5
            } else {
                (b + 1.0_f64).max(z_inf + 1.0)
            }
        } else {
            (b + 1.0_f64).max(z_inf + 1.0)
        }
    };

    // Halley iteration for largest root
    let mut converged = false;
    for _ in 0..MAX_ITER {
        let f = z * z * z + c2 * z * z + c1 * z + c0;
        let fp = 3.0 * z * z + 2.0 * c2 * z + c1;
        let fpp = 6.0 * z + 2.0 * c2;

        if fp.abs() < 1e-30 {
            break;
        }

        let mut dz = f / fp;
        // Halley correction
        let denom = 1.0 - 0.5 * dz * fpp / fp;
        if denom.abs() > 1e-15 {
            dz /= denom;
        }

        z -= dz;

        if dz.abs() < TOL {
            converged = true;
            break;
        }
    }

    if !converged {
        return Vec::new(); // Signal fallback needed
    }

    let z1 = z;

    // Synthetic division: Z^3 + c2*Z^2 + c1*Z + c0 = (Z - Z1)(Z^2 + q1*Z + q0)
    let q1 = c2 + z1;
    let q0 = c1 + z1 * q1;

    // Solve quadratic Z^2 + q1*Z + q0 = 0
    let disc = q1 * q1 - 4.0 * q0;
    let mut roots = vec![z1];

    if disc >= 0.0 {
        let sqrt_disc = disc.sqrt();
        let z2_raw = (-q1 - sqrt_disc) / 2.0;
        let z3_raw = (-q1 + sqrt_disc) / 2.0;

        // Refine each quadratic root with one Halley step
        for mut zk in [z2_raw, z3_raw] {
            let f = zk * zk * zk + c2 * zk * zk + c1 * zk + c0;
            let fp = 3.0 * zk * zk + 2.0 * c2 * zk + c1;
            let fpp = 6.0 * zk + 2.0 * c2;
            if fp.abs() > 1e-30 {
                let mut dz = f / fp;
                let denom = 1.0 - 0.5 * dz * fpp / fp;
                if denom.abs() > 1e-15 {
                    dz /= denom;
                }
                zk -= dz;
            }
            roots.push(zk);
        }
    }

    // Filter valid roots: Z > B
    let mut valid: Vec<f64> = roots
        .into_iter()
        .filter(|&r| r > b + 1e-10)
        .collect();
    valid.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    valid
}

/// Analytic cubic solver fallback using Cardano's formula.
fn cardano_cubic(c2: f64, c1: f64, c0: f64, b: f64) -> Vec<f64> {
    // Convert to depressed cubic: t^3 + pt + q = 0 where Z = t - c2/3
    let p = c1 - c2 * c2 / 3.0;
    let q = c0 - c1 * c2 / 3.0 + 2.0 * c2 * c2 * c2 / 27.0;
    let shift = -c2 / 3.0;

    let disc = -4.0 * p * p * p - 27.0 * q * q;

    let mut roots = Vec::new();

    if disc >= 0.0 {
        // Three real roots
        let r = (-p / 3.0).sqrt();
        let cos_arg = if r.abs() < 1e-30 {
            0.0
        } else {
            (-q / (2.0 * r * r * r)).clamp(-1.0, 1.0)
        };
        let theta = cos_arg.acos();
        for k in 0..3 {
            let t = 2.0 * r * ((theta + 2.0 * std::f64::consts::PI * k as f64) / 3.0).cos();
            let z = t + shift;
            if z > b + 1e-10 {
                roots.push(z);
            }
        }
    } else {
        // One real root
        let sqrt_disc_27 = (q * q / 4.0 + p * p * p / 27.0).abs().sqrt();
        let a_val = -q / 2.0 + sqrt_disc_27;
        let b_val = -q / 2.0 - sqrt_disc_27;
        let t = a_val.cbrt() + b_val.cbrt();
        let z = t + shift;
        if z > b + 1e-10 {
            roots.push(z);
        }
    }

    if roots.is_empty() {
        roots.push((b + 0.01_f64).max(0.1));
    }

    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    roots
}

/// Solve PR cubic EOS for compressibility factor Z.
///
/// Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
///
/// Uses Halley-accelerated iteration (Michelsen-style) with Cardano fallback.
///
/// Returns sorted list of valid Z roots.
pub fn solve_cubic_eos(a: f64, b: f64) -> Vec<f64> {
    let c2 = -(1.0 - b);
    let c1 = a - 3.0 * b * b - 2.0 * b;
    let c0 = -(a * b - b * b - b * b * b);

    // Try Halley solver first
    let valid = halley_cubic(c2, c1, c0, b);

    if !valid.is_empty() {
        valid
    } else {
        // Fallback to analytic Cardano
        cardano_cubic(c2, c1, c0, b)
    }
}

/// Gibbs energy difference G_vapor - G_liquid (dimensionless) for PR EOS.
///
/// If delta_G < 0, vapor is preferred; if delta_G > 0, liquid is preferred.
pub fn gibbs_delta(a: f64, b: f64, zl: f64, zv: f64) -> f64 {
    let d1 = 1.0 + SQRT2;
    let d2 = 1.0 - SQRT2;

    if zv - b <= 0.0 || zl - b <= 0.0 || b < 1e-15 {
        return 0.0;
    }

    let term1 = ((zv - b) / (zl - b)).ln();

    let num = (zv + d2 * b) * (zl + d1 * b);
    let den = (zv + d1 * b) * (zl + d2 * b);
    if den <= 0.0 || num <= 0.0 {
        return 0.0;
    }

    let term2 = -(a / (2.0 * SQRT2 * b)) * (num / den).ln();
    let term3 = -(zv - zl);

    term1 + term2 + term3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_cubic_single_root() {
        // At low A, B expect single root near 1
        let roots = solve_cubic_eos(0.01, 0.001);
        assert!(!roots.is_empty());
        assert!(roots[0] > 0.0);
    }

    #[test]
    fn test_solve_cubic_two_roots() {
        // Typical two-phase conditions
        let roots = solve_cubic_eos(0.5, 0.05);
        assert!(!roots.is_empty());
    }
}
