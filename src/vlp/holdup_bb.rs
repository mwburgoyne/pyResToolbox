/// Beggs & Brill flow pattern map, horizontal holdup, inclination correction,
/// and two-phase friction.
/// Direct port of nodal.py _bb_flow_pattern, _bb_horizontal_holdup,
/// _bb_inclination_correction and _bb_two_phase_friction.

use super::constants::*;
use super::pvt_helpers::clamp;

pub const BB_SEGREGATED: i32 = 0;
pub const BB_INTERMITTENT: i32 = 1;
pub const BB_DISTRIBUTED: i32 = 2;
pub const BB_TRANSITION: i32 = 3;

/// Flow pattern detection. Returns (pattern, transition_factor).
pub fn bb_flow_pattern(froude: f64, lambda_l: f64) -> (i32, f64) {
    if lambda_l <= 0.0 || lambda_l >= 1.0 {
        return (BB_DISTRIBUTED, 0.0);
    }
    let l1 = BB_L1_A * lambda_l.powf(BB_L1_B);
    let l2 = BB_L2_A * lambda_l.powf(BB_L2_B);
    let l3 = BB_L3_A * lambda_l.powf(BB_L3_B);
    let l4 = BB_L4_A * lambda_l.powf(BB_L4_B);
    // Revised Beggs & Brill flow pattern map (Brill and Beggs, Two-Phase Flow
    // in Pipes, 1991 edition):
    //   Segregated:   (lambda_l <  0.01 and Fr < L1) or (lambda_l >= 0.01 and Fr < L2)
    //   Transition:    lambda_l >= 0.01 and L2 <= Fr <= L3
    //   Intermittent: (0.01 <= lambda_l < 0.4 and L3 < Fr <= L1)
    //                 or (lambda_l >= 0.4 and L3 < Fr <= L4)
    //   Distributed:  (lambda_l < 0.4 and Fr >= L1) or (lambda_l >= 0.4 and Fr > L4)
    if (lambda_l < 0.01 && froude < l1) || (lambda_l >= 0.01 && froude < l2) {
        return (BB_SEGREGATED, 0.0);
    }
    if lambda_l >= 0.01 && l2 <= froude && froude <= l3 {
        let a = if l3 > l2 { (l3 - froude) / (l3 - l2) } else { 0.5 };
        return (BB_TRANSITION, clamp(a, 0.0, 1.0));
    }
    if ((0.01..0.4).contains(&lambda_l) && l3 < froude && froude <= l1)
        || (lambda_l >= 0.4 && l3 < froude && froude <= l4)
    {
        return (BB_INTERMITTENT, 0.0);
    }
    if (lambda_l < 0.4 && froude >= l1) || (lambda_l >= 0.4 && froude > l4) {
        return (BB_DISTRIBUTED, 0.0);
    }
    (BB_INTERMITTENT, 0.0)
}

/// BB horizontal holdup (no inclination correction).
pub fn bb_horizontal_holdup(lambda_l: f64, froude: f64, pattern: i32) -> f64 {
    if lambda_l <= 0.0 {
        return 0.0;
    }
    if lambda_l >= 1.0 {
        return 1.0;
    }
    if froude <= 0.0 {
        return lambda_l;
    }
    let hl0 = if pattern == BB_SEGREGATED {
        let (a, b, c) = BB_HL_SEG;
        a * lambda_l.powf(b) / froude.powf(c)
    } else if pattern == BB_INTERMITTENT {
        let (a, b, c) = BB_HL_INT;
        a * lambda_l.powf(b) / froude.powf(c)
    } else if pattern == BB_DISTRIBUTED {
        let (a, b, c) = BB_HL_DIS;
        a * lambda_l.powf(b) / froude.powf(c)
    } else {
        lambda_l
    };
    hl0.max(lambda_l)
}

/// BB inclination correction.
pub fn bb_inclination_correction(
    hl0: f64, lambda_l: f64, n_lv: f64, froude: f64,
    pattern: i32, theta: f64,
) -> f64 {
    if lambda_l <= 0.0 || lambda_l >= 1.0 {
        return hl0;
    }
    let (e_p, f_p, g_p, h_p) = if pattern == BB_SEGREGATED {
        BB_IC_SEG
    } else if pattern == BB_INTERMITTENT {
        BB_IC_INT
    } else {
        return hl0; // Distributed or Transition - no correction
    };

    let arg = e_p * lambda_l.powf(f_p)
        * n_lv.max(1e-10).powf(g_p)
        * froude.max(1e-10).powf(h_p);
    let mut c_corr = if arg > 0.0 {
        (1.0 - lambda_l) * arg.ln()
    } else {
        0.0
    };
    c_corr = c_corr.max(0.0);
    let sin18 = (1.8 * theta).sin();
    let psi = 1.0 + c_corr * (sin18 - 0.333 * sin18.powi(3));
    clamp(hl0 * psi, lambda_l, 1.0)
}

/// BB two-phase friction multiplier.
pub fn bb_two_phase_friction(f_ns: f64, lambda_l: f64, hl_theta: f64) -> f64 {
    if hl_theta <= 0.0 {
        return f_ns;
    }
    let y = lambda_l / (hl_theta * hl_theta);
    let s = if y <= 0.0 || y == 1.0 {
        0.0
    } else if 1.0 < y && y < 1.2 {
        (2.2 * y - 1.2).ln()
    } else {
        let ln_y = y.ln();
        let denom = BB_SF_C0 + BB_SF_C1 * ln_y + BB_SF_C2 * ln_y * ln_y
            + BB_SF_C3 * ln_y.powi(4);
        if denom.abs() >= 1e-6 { ln_y / denom } else { 0.0 }
    };
    let s = clamp(s, BB_SF_LO, BB_SF_HI);
    f_ns * s.exp()
}
