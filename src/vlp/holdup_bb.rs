/// Beggs & Brill flow pattern map, horizontal holdup, inclination correction, and two-phase friction.
/// Direct port of nodal.py lines 1558-1642.

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
    let l1 = 316.0 * lambda_l.powf(0.302);
    let l2 = 0.0009252 * lambda_l.powf(-2.4684);
    let l3 = 0.10 * lambda_l.powf(-1.4516);
    let l4 = 0.5 * lambda_l.powf(-6.738);

    if lambda_l < 0.01 && froude < l1 {
        return (BB_SEGREGATED, 0.0);
    }
    if lambda_l >= 0.01 && froude < l2 {
        return (BB_SEGREGATED, 0.0);
    }
    if lambda_l >= 0.01 && l2 <= froude && froude <= l3 {
        let a = if l3 > l2 { (l3 - froude) / (l3 - l2) } else { 0.5 };
        return (BB_TRANSITION, clamp(a, 0.0, 1.0));
    }
    if (lambda_l >= 0.01 && froude > l3 && froude < l1)
        || (lambda_l < 0.01 && froude >= l1)
    {
        return (BB_INTERMITTENT, 0.0);
    }
    if lambda_l < 0.4 && froude >= l1 {
        return (BB_DISTRIBUTED, 0.0);
    }
    if lambda_l >= 0.4 && froude > l4 {
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
        0.98 * lambda_l.powf(0.4846) / froude.powf(0.0868)
    } else if pattern == BB_INTERMITTENT {
        0.845 * lambda_l.powf(0.5351) / froude.powf(0.0173)
    } else if pattern == BB_DISTRIBUTED {
        1.065 * lambda_l.powf(0.5824) / froude.powf(0.0609)
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
        (0.011, -3.7680, 3.5390, -1.6140)
    } else if pattern == BB_INTERMITTENT {
        (2.960, 0.3050, -0.4473, 0.0978)
    } else {
        return hl0; // Distributed or Transition — no correction
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
        let denom = -0.0523 + 3.182 * ln_y - 0.8725 * ln_y * ln_y
            + 0.01853 * ln_y.powi(4);
        if denom.abs() >= 1e-6 { ln_y / denom } else { 0.0 }
    };
    let s = clamp(s, -5.0, 5.0);
    f_ns * s.exp()
}
