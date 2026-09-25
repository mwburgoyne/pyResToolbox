//! Salt legs of the gas-free brine viscosity chain, specialised to a pure NaCl
//! brine of molality m (mol/kg water):
//!
//! * the ion-additive modified Jones-Dole ratio of PHREEQC (Appelo &
//!   Parkhurst), port of brine/jones_dole_viscosity.py viscosity_ratio;
//! * Kestin, Khalifa & Correia (1981) JPCRD 10, 71 NaCl viscosity, used only
//!   for its pressure factor, port of brine/kestin_nacl_viscosity.py and
//!   brine/viscosity_route.py pressure_factor.
//!
//! For NaCl the general multi-ion code collapses: I = m, z1 = z2 = 1, Cl- is
//! the Jones-Dole reference solute (B = D = 0), and the anion-volume factor
//! fan = 2 - V_Cl/V_Cl = 1, so the Vm_ion chain is not needed.

use super::water::dielectric_constant;

// Parameters from PHREEQC pitzer.dat as embedded in brine/_lib_pitzer_params.py.
/// Na+ `-viscosity` line: b0, b1, b2, d1, d2, d3 (slots 7, 8 are zero).
const JD_NA: [f64; 6] = [0.1387, -0.0866, 0.0125, 0.0145, 0.0075, 1.062];
/// Tracer diffusion coefficients (m2/s) and their temperature parameters.
const DW_NA: (f64, f64) = (1.33e-9, 75.0);
const DW_CL: (f64, f64) = (2.033e-9, 216.0);

/// PHREEQC transport.cpp: pure-water viscosity at 25 degC, mPa s.
const VISCOS_0_25: f64 = 0.8900239182946;
/// transport.cpp caps tc for the viscosity model.
const TC_MAX: f64 = 200.0;
/// Falkenhagen-Dole prefactor, transport.cpp.
const FD_PREFACTOR: f64 = 4.3787e-14;

/// B*m and D*m terms for Na+ (z = 1). Port of jones_dole_viscosity._B_and_D.
fn na_b_and_d(tc: f64, m: f64, ionic: f64) -> (f64, f64) {
    let [b0, b1, b2, d1, d2, d3] = JD_NA;
    let b_term = (b0 + b1 * (-b2 * tc).exp()) * m;
    let f_z = 1.0; // (z^2 + |z|)/2 for z = 1
    let f_i = if d3 >= 1.0 {
        ionic / 3.0 / d3
    } else if d3 > 0.4 {
        -0.8 / d3
    } else {
        -1.0
    };
    let mut d_term = ((d1 * (-d2 * tc).exp()) * m
        * (ionic.powf(d3) * (1.0 + f_i) + (m * f_z).powf(d3)))
        / (2.0 + f_i);
    if d_term < -1e-5 {
        d_term = 0.0;
    }
    (b_term, d_term)
}

/// Tracer diffusion coefficient scaled to mu_0 and T, transport.cpp.
fn scaled_dw(dw: (f64, f64), t_k: f64, scale: f64) -> f64 {
    let mut d = dw.0 * scale;
    if dw.1 != 0.0 {
        d *= (dw.1 / t_k - dw.1 / 298.15).exp();
    }
    d
}

/// Falkenhagen-Dole electrostatic term times sqrt(eq/2), already multiplied
/// by mu_0. Port of jones_dole_viscosity._falkenhagen_dole_A for NaCl.
fn falkenhagen_dole_a(t_k: f64, p_mpa: f64, m: f64, mu_0: f64) -> f64 {
    let eps_r = dielectric_constant(t_k, p_mpa * 10.0);
    let scale = VISCOS_0_25 / mu_0;
    let dw_plus = scaled_dw(DW_NA, t_k, scale);
    let dw_min = scaled_dw(DW_CL, t_k, scale);
    // Accumulators exactly as the Python loop builds them (z = +1 and -1).
    let (m_plus, eq_plus, eq_dw_plus) = (m, m, m / dw_plus);
    let (m_min, eq_min, eq_dw_min) = (m, m, m / dw_min);
    if m_plus == 0.0 || m_min == 0.0 || eq_dw_plus == 0.0 || eq_dw_min == 0.0 {
        return 0.0;
    }
    let z1 = eq_plus / m_plus;
    let z2 = eq_min / m_min;
    let d1 = eq_plus / eq_dw_plus;
    let d2 = eq_min / eq_dw_min;
    let t1 = (d1 - d2) / ((d1 * z1 + d2 * z2).sqrt() + ((d1 + d2) * (z1 + z2)).sqrt());
    let psi = (d1 * z2 + d2 * z1) / 4.0 - z1 * z2 * t1 * t1;
    let a = FD_PREFACTOR * t_k.powf(1.5)
        / ((eps_r * (z1 + z2) / z1.max(z2)).sqrt() * (d1 * d2))
        * psi;
    a * ((eq_plus + eq_min) / 2.0).sqrt()
}

/// Jones-Dole mu(brine)/mu(water) for NaCl molality m > 0, given mu_0 (mPa s)
/// at the same T and P. Port of jones_dole_viscosity.viscosity_ratio.
pub fn jones_dole_ratio(t_k: f64, p_mpa: f64, m: f64, mu_0: f64) -> f64 {
    let tc = (t_k - 273.15).min(TC_MAX);
    let ionic = m; // 0.5 * (m * 1 + m * 1)
    let (bc, mut dc) = if m <= 1e-9 { (0.0, 0.0) } else { na_b_and_d(tc, m, ionic) };
    if dc < 0.0 {
        dc = 0.0;
    }
    let a_term = falkenhagen_dole_a(t_k, p_mpa, m, mu_0);
    let fan = 1.0;
    1.0 + a_term / mu_0 + fan * (bc + dc)
}

// ----------------------------------------------------------------------------
// Kestin, Khalifa & Correia (1981), Eqs. (1)-(10). t in degC, p in MPa.
// ----------------------------------------------------------------------------

const K_ALPHA: [f64; 4] = [1.2378, -1.303e-3, 3.06e-6, 2.55e-8]; // Eq. (3)
const K_MU_W0_20C: f64 = 1002.0; // micro Pa s
const K_A: [f64; 3] = [3.324e-2, 3.624e-3, -1.879e-4]; // Eq. (4)
const K_B: [f64; 3] = [-3.96e-2, 1.02e-2, -7.02e-4]; // Eq. (5)
const K_BETA_W: [f64; 5] = [-1.297, 5.74e-2, -6.97e-4, 4.47e-6, -1.05e-8]; // Eq. (7)
const K_GAMMA: [f64; 2] = [0.545, 2.8e-3]; // Eq. (8)
const K_MS: [f64; 3] = [6.044, 2.8e-3, 3.6e-5]; // Eq. (9)
const K_BETA_STAR: [f64; 3] = [2.5, -2.0, 0.5]; // Eq. (10)

/// Kestin's calibration box; the pressure factor is clamped to it
/// (viscosity_route._P_FACTOR_*).
const PF_T_C: (f64, f64) = (20.0, 150.0);
const PF_P_MPA: (f64, f64) = (0.1, 35.0);
const PF_M_MAX: f64 = 6.0;

/// sum_i c_i * x^(i + offset), accumulated in index order as Python's sum().
fn poly(coef: &[f64], x: f64, offset: i32) -> f64 {
    let mut s = 0.0;
    for (i, c) in coef.iter().enumerate() {
        s += c * x.powi(i as i32 + offset);
    }
    s
}

fn kestin_mu(t_c: f64, p_mpa: f64, m: f64) -> f64 {
    let muw = K_MU_W0_20C * 10.0_f64.powf(poly(&K_ALPHA, 20.0 - t_c, 1) / (96.0 + t_c));
    let mu0 = muw * 10.0_f64.powf(poly(&K_A, m, 1) + poly(&K_B, m, 1) * (muw / K_MU_W0_20C).log10());
    let bw = poly(&K_BETA_W, t_c, 0);
    let beta_es = K_GAMMA[0] + K_GAMMA[1] * t_c - bw;
    let r = m / poly(&K_MS, t_c, 0);
    let beta = beta_es * poly(&K_BETA_STAR, r, 1) + bw;
    mu0 * (1.0 + beta * p_mpa / 1000.0)
}

fn kestin_salt_ratio(t_c: f64, p_mpa: f64, m: f64) -> f64 {
    kestin_mu(t_c, p_mpa, m) / kestin_mu(t_c, p_mpa, 0.0)
}

/// Kestin's pressure dependence of the salt ratio normalised to 0.1 MPa, at
/// ionic strength m, clamped to the calibration box. Port of
/// viscosity_route.pressure_factor for NaCl (m > 0).
pub fn kestin_pressure_factor(t_k: f64, p_mpa: f64, m: f64) -> f64 {
    let t = (t_k - 273.15).clamp(PF_T_C.0, PF_T_C.1);
    let p = p_mpa.clamp(PF_P_MPA.0, PF_P_MPA.1);
    let mm = m.min(PF_M_MAX);
    kestin_salt_ratio(t, p, mm) / kestin_salt_ratio(t, 0.1, mm)
}
