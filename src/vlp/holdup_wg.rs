/// Woldesemayat-Ghajar drift-flux void fraction and LM friction gradient.
/// Direct port of nodal.py lines 1044-1110.

use super::pvt_helpers::clamp;
use super::friction::serghides_fanning;

const G_SI: f64 = 9.80665;
const P_ATM_PA: f64 = 101325.0;

/// WG drift-flux void fraction (gas fraction).
pub fn wg_void_fraction(
    u_sg: f64, u_sl: f64, rho_g: f64, rho_l: f64,
    sigma: f64, diam: f64, theta: f64, p_sys: f64,
) -> f64 {
    if u_sg <= 0.0 {
        return 0.0;
    }
    if u_sl <= 0.0 {
        return 1.0;
    }
    if rho_l <= rho_g || sigma <= 0.0 || diam <= 0.0 || p_sys <= 0.0 {
        return 0.5;
    }
    let u_m = u_sg + u_sl;
    let lam = u_sg / u_m;
    let density_ratio = (rho_g / rho_l).powf(0.1);
    let co = lam * (1.0 + (u_sl / u_sg).powf(density_ratio));
    let buoy = G_SI * diam * sigma * (1.0 + theta.cos())
        * (rho_l - rho_g) / (rho_l * rho_l);
    let u_gm_base = 2.9 * buoy.max(0.0).powf(0.25);
    let p_exp = P_ATM_PA / p_sys;
    let incl_factor = (1.22 + 1.22 * theta.sin()).powf(p_exp);
    let u_gm = u_gm_base * incl_factor;
    let denom = co * u_m + u_gm;
    if denom <= 0.0 {
        return 0.0;
    }
    clamp(u_sg / denom, 0.0, 1.0)
}

/// Lockhart-Martinelli two-phase friction gradient (Pa/m).
pub fn wg_friction_gradient_lm(
    m_flow_g: f64, m_flow_l: f64, rho_g: f64, rho_l: f64,
    mu_g: f64, mu_l: f64, diam: f64, rough: f64,
) -> f64 {
    if diam <= 0.0 {
        return 0.0;
    }
    let area_si = std::f64::consts::PI * diam * diam / 4.0;
    if area_si <= 0.0 {
        return 0.0;
    }
    let m_total = m_flow_g + m_flow_l;
    if m_total < 1e-30 {
        return 0.0;
    }
    let mass_flux = m_total / area_si;
    let x = m_flow_g / m_total;
    let e_over_d = rough / diam;

    if x <= 0.0 {
        let re_l = mass_flux * diam / mu_l;
        let f_l = serghides_fanning(re_l, e_over_d);
        return 2.0 * f_l * mass_flux * mass_flux / (diam * rho_l);
    }
    if x >= 1.0 {
        let re_g = mass_flux * diam / mu_g;
        let f_g = serghides_fanning(re_g, e_over_d);
        return 2.0 * f_g * mass_flux * mass_flux / (diam * rho_g);
    }

    let g_l = mass_flux * (1.0 - x);
    let g_g = mass_flux * x;
    let re_l = g_l * diam / mu_l;
    let re_g = g_g * diam / mu_g;
    let f_l = serghides_fanning(re_l, e_over_d);
    let f_g = serghides_fanning(re_g, e_over_d);
    let dpdz_l = 2.0 * f_l * g_l * g_l / (diam * rho_l);
    let dpdz_g = 2.0 * f_g * g_g * g_g / (diam * rho_g.max(1e-30));
    let x_param = (dpdz_l / dpdz_g.max(1e-30)).sqrt();
    if x_param < 1e-30 {
        return dpdz_g;
    }
    let liq_turb = re_l >= 2100.0;
    let gas_turb = re_g >= 2100.0;
    let chisholm_c = match (liq_turb, gas_turb) {
        (true, true) => 20.0,
        (false, true) => 12.0,
        (true, false) => 10.0,
        (false, false) => 5.0,
    };
    let phi2_l = 1.0 + chisholm_c / x_param + 1.0 / (x_param * x_param);
    phi2_l * dpdz_l
}
