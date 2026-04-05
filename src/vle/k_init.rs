/// K-value initialization for VLE flash.
///
/// S&W-specific K-value correlations fitted to binary VLE results,
/// replacing the standard Wilson correlation which gives K < 1 for
/// many gas-water systems at moderate conditions.

use crate::vle::components::*;

// =============================================================================
// S&W K-value parameters
// =============================================================================

/// Light gases: Cross form [a, b, c, d, e, f]
/// ln(K) = a + b(Tc/T) + c·ln(P/Pc) + d(Tc/T)^2 + e·ln(P/Pc)^2 + f(Tc/T)·ln(P/Pc)
const SW_KVALUE_PARAMS: [(usize, [f64; 6]); 9] = [
    (IDX_H2,     [6.4295, 25.5844, -0.5985, -50.0000, -0.0007, -3.2598]),
    (IDX_CO2,    [-5.9974, 26.3804, -0.9380, -16.3941, 0.0607, 0.3688]),
    (IDX_N2,     [1.4998, 36.4929, -0.6419, -50.0000, 0.0110, -0.6328]),
    (IDX_H2S,    [-1.0549, 7.7334, -1.3646, -3.4035, 0.0850, 0.8651]),
    (IDX_CH4,    [-2.7107, 36.5276, -0.7040, -33.2023, 0.0312, -0.1419]),
    (IDX_C2H6,   [-9.8804, 40.0208, -0.9108, -22.8267, 0.0756, 0.4366]),
    (IDX_C3H8,   [-9.8805, 33.3750, -1.0803, -15.1332, 0.0836, 0.6959]),
    (IDX_IC4H10, [-8.8362, 29.1657, -1.0755, -11.0413, 0.0685, 0.7280]),
    (IDX_NC4H10, [-8.4159, 27.0161, -1.0886, -10.2525, 0.0613, 0.7515]),
];

/// Heavy HCs: LogLinear form [a, b, c] + floor=10
const SW_KVALUE_HEAVY: [(usize, [f64; 3]); 6] = [
    (IDX_IC5H12,  [6.5325, 3.4733, -0.0386]),
    (IDX_NC5H12,  [6.2726, 3.7905, -0.0083]),
    (IDX_NC6H14,  [4.5498, 6.4018, 0.1051]),
    (IDX_NC7H16,  [2.6475, 8.9734, 0.1954]),
    (IDX_NC8H18,  [0.3846, 11.9724, 0.2711]),
    (IDX_NC10H22, [3.8490, -1.1163, -0.0764]),
];

/// Water K-value: ln(K_H2O) = a + b/T + c·ln(P) + d/T^2 + e·ln(P)/T + f·ln(P)^2
const SW_KVALUE_WATER: [f64; 6] = [
    34.3692, -1873.3949, -3.2671, -555202.1539, 42.1759, 0.0804,
];

/// Find cross-form parameters for a light gas component index.
fn find_light_params(comp_idx: usize) -> Option<[f64; 6]> {
    for &(idx, params) in &SW_KVALUE_PARAMS {
        if idx == comp_idx {
            return Some(params);
        }
    }
    None
}

/// Find heavy HC parameters for a component index.
fn find_heavy_params(comp_idx: usize) -> Option<[f64; 3]> {
    for &(idx, params) in &SW_KVALUE_HEAVY {
        if idx == comp_idx {
            return Some(params);
        }
    }
    None
}

/// S&W-specific K-value initialization for gas-water flash.
///
/// Returns array of K-values for all components in the given index order.
///
/// # Arguments
/// * `comp_indices` - Component index array
/// * `tc` - Critical temperatures (K)
/// * `pc` - Critical pressures (Pa)
/// * `omega` - Acentric factors
/// * `t_k` - Temperature (K)
/// * `p_pa` - Pressure (Pa)
pub fn sw_kvalue_init(
    comp_indices: &[usize],
    tc: &[f64],
    pc: &[f64],
    omega: &[f64],
    t_k: f64,
    p_pa: f64,
) -> Vec<f64> {
    let nc = comp_indices.len();
    let mut k_vals = vec![1.0; nc];
    let t = t_k;
    let p = p_pa;

    for i in 0..nc {
        let idx = comp_indices[i];

        if idx == IDX_H2O {
            // Universal water K-value (6-param T-P fit)
            let a = SW_KVALUE_WATER;
            let ln_p = p.ln();
            k_vals[i] = (a[0] + a[1] / t + a[2] * ln_p + a[3] / (t * t)
                + a[4] * ln_p / t + a[5] * ln_p * ln_p)
                .exp();
        } else if let Some(params) = find_light_params(idx) {
            // Light gas: Cross form
            let tr_inv = tc[i] / t;
            let ln_pr = (p / pc[i]).ln();
            k_vals[i] = (params[0]
                + params[1] * tr_inv
                + params[2] * ln_pr
                + params[3] * tr_inv * tr_inv
                + params[4] * ln_pr * ln_pr
                + params[5] * tr_inv * ln_pr)
                .exp();
        } else if let Some(params) = find_heavy_params(idx) {
            // Heavy HC: LogLinear + floor
            let tr_inv = tc[i] / t;
            let ln_pr = (p / pc[i]).ln();
            k_vals[i] = (params[0] + params[1] * tr_inv + params[2] * ln_pr)
                .exp()
                .max(10.0);
        } else {
            // Fallback: standard Wilson for unknown components
            k_vals[i] = (pc[i] / p)
                * (5.373 * (1.0 + omega[i]) * (1.0 - tc[i] / t)).exp();
        }
    }

    // Clip to [1e-10, 1e10]
    for v in k_vals.iter_mut() {
        *v = v.clamp(1e-10, 1e10);
    }

    k_vals
}

/// Standard Wilson K-value correlation.
///
/// K_i = (Pc_i / P) * exp(5.373 * (1 + omega_i) * (1 - Tc_i / T))
#[allow(dead_code)]
pub fn wilson_k_value(tc: f64, pc: f64, omega: f64, t_k: f64, p_pa: f64) -> f64 {
    (pc / p_pa) * (5.373 * (1.0 + omega) * (1.0 - tc / t_k)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sw_kvalue_water() {
        // At moderate conditions, water K should be << 1
        let comp = [IDX_H2O, IDX_CH4];
        let tc = [647.3, 190.6];
        let pc = [22.12e6, 4.60e6];
        let omega = [0.3434, 0.0108];
        let k = sw_kvalue_init(&comp, &tc, &pc, &omega, 373.15, 100e5);
        assert!(k[0] < 1.0, "Water K should be < 1, got {}", k[0]);
        assert!(k[1] > 1.0, "CH4 K should be > 1, got {}", k[1]);
    }

    #[test]
    fn test_sw_kvalue_clamp() {
        // All K values should be within bounds
        let comp = [IDX_H2O, IDX_H2];
        let tc = [647.3, 33.145];
        let pc = [22.12e6, 1.2964e6];
        let omega = [0.3434, -0.219];
        let k = sw_kvalue_init(&comp, &tc, &pc, &omega, 300.0, 1e5);
        for &ki in &k {
            assert!(ki >= 1e-10 && ki <= 1e10);
        }
    }
}
