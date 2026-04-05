/// Simplified PVT helpers for VLP segment loops.
/// Direct port of nodal.py lines 430-582.

const MW_AIR: f64 = 28.97;
const LN10: f64 = std::f64::consts::LN_10;

#[inline]
fn log10_safe(x: f64) -> f64 {
    if x <= 0.0 { -30.0 } else { x.ln() / LN10 }
}

#[inline]
pub fn clamp(val: f64, lo: f64, hi: f64) -> f64 {
    if val < lo { lo } else if val > hi { hi } else { val }
}

/// Sutton critical properties from gas SG.
#[inline]
pub fn sutton_tc_pc(sg: f64) -> (f64, f64) {
    let tpc = 169.2 + 349.5 * sg - 74.0 * sg * sg;
    let ppc = 756.8 - 131.0 * sg - 3.6 * sg * sg;
    (tpc, ppc)
}

/// Z-factor via Hall-Yarborough (1973).
pub fn z_factor(_sg: f64, temp_f: f64, press_psia: f64, tc: f64, pc: f64) -> f64 {
    if press_psia < 1.0 {
        return 1.0;
    }
    let tr = (temp_f + 459.67) / tc;
    let pr = press_psia / pc;

    let t_inv = 1.0 / tr;
    let a = -0.06125 * t_inv * (-1.2 * (1.0 - t_inv).powi(2)).exp();
    let b = 14.76 * t_inv - 9.76 * t_inv.powi(2) + 4.58 * t_inv.powi(3);
    let c = 90.7 * t_inv - 242.2 * t_inv.powi(2) + 42.4 * t_inv.powi(3);
    let d = 2.18 + 2.82 * t_inv;

    let mut y = clamp(0.0125 * pr * t_inv, 1e-10, 0.9);

    for _ in 0..50 {
        let y2 = y * y;
        let y3 = y2 * y;
        let y4 = y2 * y2;
        let one_m_y = 1.0 - y;
        if one_m_y.abs() < 1e-15 {
            break;
        }
        let fy = (y + y2 + y3 - y4) / one_m_y.powi(3)
            + a * pr - b * y2 + c * y.powf(d);
        let dfy = (1.0 + 4.0 * y + 4.0 * y2 - 4.0 * y3 + y4) / one_m_y.powi(4)
            - 2.0 * b * y + c * d * y.powf(d - 1.0);
        if dfy.abs() < 1e-30 {
            break;
        }
        let dy = fy / dfy;
        y -= dy;
        y = clamp(y, 1e-10, 0.99);
        if dy.abs() < 1e-12 {
            break;
        }
    }

    if y.abs() < 1e-30 {
        return 1.0;
    }
    let z = -a * pr / y;
    clamp(z, 0.1, 5.0)
}

/// Gas viscosity via Lee-Gonzalez-Eakin (1966) (cP).
pub fn gas_viscosity(sg: f64, temp_f: f64, press_psia: f64, z: f64, _tc: f64, _pc: f64) -> f64 {
    let temp_r = temp_f + 459.67;
    let mw = MW_AIR * sg;
    let rho_gcc = (MW_AIR * sg * press_psia / (z * 10.732 * temp_r)) / 62.428;

    let k_val = (9.4 + 0.02 * mw) * temp_r.powf(1.5) / (209.0 + 19.0 * mw + temp_r);
    let x_val = 3.5 + 986.0 / temp_r + 0.01 * mw;
    let y_val = 2.4 - 0.2 * x_val;
    let mut exp_arg = x_val * rho_gcc.powf(y_val);
    if exp_arg > 500.0 {
        exp_arg = 500.0;
    }
    (1e-4 * k_val * exp_arg.exp()).max(1e-6)
}

/// Gas density in lb/ft^3.
#[inline]
#[allow(dead_code)]
pub fn gas_density(sg: f64, temp_f: f64, press_psia: f64, z: f64) -> f64 {
    MW_AIR * sg * press_psia / (z * 10.732 * (temp_f + 459.67))
}

/// Simplified water viscosity (cP).
pub fn water_viscosity(press_psia: f64, temp_f: f64, _salinity: f64) -> f64 {
    let t = if temp_f < 32.0 { 32.0 } else { temp_f };
    let a_coeff = -3.79418 + 604.129 / (139.18 + t);
    let mut mu_w = 10.0_f64.powf(a_coeff);
    // salinity correction omitted (VLP helpers pass 0)
    mu_w *= 1.0 + 5e-4 * (press_psia - 14.7) / 1000.0;
    mu_w.max(0.01)
}

/// Standing (1947) solution GOR estimate (scf/STB).
pub fn standing_rs(gsg: f64, press_psia: f64, temp_f: f64, api: f64) -> f64 {
    let exponent = 0.00091 * temp_f - 0.0125 * api;
    let denom = 18.2 * 10.0_f64.powf(exponent);
    if denom <= 0.0 {
        return 0.0;
    }
    gsg * (press_psia / denom).max(0.0).powf(1.2048)
}

/// Velarde (1997) solution GOR vs pressure (scf/STB).
pub fn velarde_rs(sgsp: f64, api: f64, temp_f: f64, pb: f64, rsb: f64, press_psia: f64) -> f64 {
    if press_psia >= pb {
        return rsb;
    }
    let pb_diff = (pb - 14.7_f64).max(0.1);
    let a1 = 0.000000973 * sgsp.powf(1.672608) * api.powf(0.92987)
        * temp_f.powf(0.247235) * pb_diff.powf(1.056052);
    let a2 = 0.022339 * sgsp.powf(-1.004750) * api.powf(0.337711)
        * temp_f.powf(0.132795) * pb_diff.powf(0.302065);
    let a3 = 0.725167 * sgsp.powf(-1.485480) * api.powf(-0.164741)
        * temp_f.powf(-0.09133) * pb_diff.powf(0.047094);
    let prr = ((press_psia - 14.7) / pb_diff).max(0.0);
    (a1 * prr.powf(a2) + (1.0 - a1) * prr.powf(a3)) * rsb
}

/// Reservoir oil density via McCain iterative method (lb/ft^3).
pub fn oil_density_mccain(rs: f64, sgsp: f64, sgsto: f64, press_psia: f64, temp_f: f64) -> f64 {
    let mut rhopo1 = 52.8 - 0.01 * rs;
    let mut rhopo2 = rhopo1;
    for _ in 0..100 {
        let rhoa2 = -49.893 + 85.0149 * sgsp - 3.70373 * sgsp * rhopo1
            + 0.0479818 * sgsp * rhopo1 * rhopo1 + 2.98914 * rhopo1
            - 0.035688 * rhopo1 * rhopo1;
        let denom = if rhoa2.abs() > 1e-30 {
            73.71 + rs * sgsp / rhoa2
        } else {
            73.71
        };
        rhopo2 = (rs * sgsp + 4600.0 * sgsto) / denom;
        if (rhopo1 - rhopo2).abs() < 0.0003 {
            break;
        }
        rhopo1 = rhopo2;
    }

    let rhopo = rhopo2;
    let drhop = (0.167 + 16.181 * 10.0_f64.powf(-0.0425 * rhopo))
        * (press_psia / 1000.0)
        - 0.01 * (0.299 + 263.0 * 10.0_f64.powf(-0.0603 * rhopo))
            * (press_psia / 1000.0).powi(2);
    let rhobs = rhopo + drhop;
    let dt = (temp_f - 60.0).max(0.001);
    let drhot = (0.00302 + 1.505 * rhobs.powf(-0.951)) * dt.powf(0.938)
        - (0.0216 - 0.0233 * 10.0_f64.powf(-0.0161 * rhobs)) * dt.powf(0.475);
    rhobs - drhot
}

/// Oil viscosity (Beggs-Robinson + undersaturated correction) (cP).
pub fn oil_viscosity_full(
    sgsp: f64, api: f64, temp_f: f64, rsb: f64, pb: f64,
    press_psia: f64, vis_frac: f64, rsb_frac: f64,
) -> f64 {
    let rs = velarde_rs(sgsp, api, temp_f, pb, rsb / rsb_frac, press_psia) * rsb_frac;
    let c = 10.0_f64.powf(3.0324 - 0.02023 * api) * temp_f.powf(-1.163);
    let mu_od = 10.0_f64.powf(c) - 1.0;
    let a = 10.715 * (rs + 100.0).powf(-0.515);
    let b = 5.44 * (rs + 150.0).powf(-0.338);
    let mut mu_or = a * mu_od.powf(b);

    if press_psia > pb {
        let ab = 10.715 * (rsb + 100.0).powf(-0.515);
        let bb = 5.44 * (rsb + 150.0).powf(-0.338);
        let mu_orb = ab * mu_od.powf(bb);
        let log_mu = log10_safe(mu_orb);
        let aa = -1.0146 + 1.3322 * log_mu - 0.4876 * log_mu.powi(2)
            - 1.15036 * log_mu.powi(3);
        mu_or = mu_orb + 0.0013449 * (press_psia - pb) * 10.0_f64.powf(aa);
    }

    mu_or.max(0.001) * vis_frac
}
