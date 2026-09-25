//! Pure-water legs of the gas-free brine viscosity chain: IAPWS-IF97 Region 1
//! density, IAPWS-2008 viscosity (industrial form, mu_2 = 1) and the Bradley &
//! Pitzer (1979) dielectric constant.
//!
//! Direct ports of plyasunov/iapws_if97.py rho_if97, brine/iapws_viscosity.py
//! mu_iapws2008 and brine/_lib_dielectric.py dielectric_constant. The Python
//! modules are the reference; coefficients are copied verbatim with their
//! sources. No range checks here: the caller clamps into Region 1.

// ----------------------------------------------------------------------------
// IAPWS-IF97 Region 1 (Wagner et al. 2000, Table 2). T in K, P in MPa.
// ----------------------------------------------------------------------------

const IF97_R: f64 = 461.526e-6; // MPa m3/(kg K)
const IF97_P_STAR: f64 = 16.53; // MPa
const IF97_T_STAR: f64 = 1386.0; // K

/// Region 1 (I_i, J_i, n_i). Rows with I = 0 do not enter gamma_pi and are
/// kept only so the table matches the Python source line for line.
const IF97_IJN: [(i32, i32, f64); 34] = [
    (0, -2, 0.14632971213167e+00),
    (0, -1, -0.84548187389013e+00),
    (0, 0, -0.37563603672040e+01),
    (0, 1, 0.33855169168385e+01),
    (0, 2, -0.95791963387872e+00),
    (0, 3, 0.15772038513228e+00),
    (0, 4, -0.16616417199501e-01),
    (0, 5, 0.81214629983568e-03),
    (1, -9, 0.28319080123804e-03),
    (1, -7, -0.60706301565874e-03),
    (1, -1, -0.18990068218419e-01),
    (1, 0, -0.32529748770505e-01),
    (1, 1, -0.21841717175414e-01),
    (1, 3, -0.52838357969930e-04),
    (2, -3, -0.47184321073267e-03),
    (2, 0, -0.30001780793026e-03),
    (2, 1, 0.47661393906987e-04),
    (2, 3, -0.44141845330846e-05),
    (2, 17, -0.72694996297594e-15),
    (3, -4, -0.31679644845054e-04),
    (3, 0, -0.28270797985312e-05),
    (3, 6, -0.85205128120103e-09),
    (4, -5, -0.22425281908000e-05),
    (4, -2, -0.65171222895601e-06),
    (4, 10, -0.14341729937924e-12),
    (5, -8, -0.40516996860117e-06),
    (8, -11, -0.12734301741682e-08),
    (8, -6, -0.17424871230634e-09),
    (21, -29, -0.68762131295531e-18),
    (23, -31, 0.14478307828521e-19),
    (29, -38, 0.26335781662795e-22),
    (30, -39, -0.11947622640071e-22),
    (31, -40, 0.18228094581404e-23),
    (32, -41, -0.93537087292458e-25),
];

/// Pure-water density (kg/m3), IAPWS-IF97 Region 1. Port of rho_if97.
pub fn rho_if97(t_k: f64, p_mpa: f64) -> f64 {
    let pi = p_mpa / IF97_P_STAR;
    let tau = IF97_T_STAR / t_k;
    let a = 7.1 - pi;
    let b = tau - 1.222;
    let mut gp = 0.0;
    for &(i, j, n) in IF97_IJN.iter() {
        let b_j = b.powi(j);
        if i == 0 {
            continue;
        } else if i == 1 {
            gp += -n * b_j;
        } else {
            gp += n * (-(i as f64)) * a.powi(i - 1) * b_j;
        }
    }
    IF97_P_STAR / (IF97_R * t_k * gp)
}

// ----------------------------------------------------------------------------
// IAPWS-2008 viscosity, Huber et al. (2009) JPCRD 38, 101, Eq. (36).
// ----------------------------------------------------------------------------

const VISC_T_STAR: f64 = 647.096; // K
const VISC_RHO_STAR: f64 = 322.0; // kg/m3
const VISC_MU_STAR: f64 = 1.0e-6; // Pa s

/// Huber Table 2, H_i for mu_0.
const VISC_H: [f64; 4] = [1.677_52, 2.204_62, 0.636_656_4, -0.241_605];

/// Huber Table 3, (i, j, H_ij) for mu_1, in the paper's row order.
const VISC_HIJ: [(i32, i32, f64); 21] = [
    (0, 0, 5.200_94e-1),
    (1, 0, 8.508_95e-2),
    (2, 0, -1.083_74),
    (3, 0, -2.895_55e-1),
    (0, 1, 2.225_31e-1),
    (1, 1, 9.991_15e-1),
    (2, 1, 1.887_97),
    (3, 1, 1.266_13),
    (5, 1, 1.205_73e-1),
    (0, 2, -2.813_78e-1),
    (1, 2, -9.068_51e-1),
    (2, 2, -7.724_79e-1),
    (3, 2, -4.898_37e-1),
    (4, 2, -2.570_40e-1),
    (0, 3, 1.619_13e-1),
    (1, 3, 2.573_99e-1),
    (0, 4, -3.253_72e-2),
    (3, 4, 6.984_52e-2),
    (4, 5, 8.721_02e-3),
    (3, 6, -4.356_73e-3),
    (5, 6, -5.932_64e-4),
];

/// Water viscosity (Pa s) from T (K) and density (kg/m3). Port of mu_iapws2008.
pub fn mu_iapws2008(t_k: f64, rho: f64) -> f64 {
    let tbar = t_k / VISC_T_STAR;
    let mut denom = 0.0;
    let mut tpow = 1.0;
    for &h in VISC_H.iter() {
        denom += h / tpow;
        tpow *= tbar;
    }
    let mu0 = 100.0 * tbar.sqrt() / denom;

    let rhobar = rho / VISC_RHO_STAR;
    let d_t = 1.0 / tbar - 1.0;
    let d_r = rhobar - 1.0;
    let mut total = 0.0;
    for &(i, j, h) in VISC_HIJ.iter() {
        total += h * d_t.powi(i) * d_r.powi(j);
    }
    let mu1 = (rhobar * total).exp();
    mu0 * mu1 * VISC_MU_STAR
}

/// Pure-water viscosity in mPa s (cP), IAPWS-2008 with IF97 density. Port of
/// jones_dole_viscosity.mu_water.
pub fn mu_water(t_k: f64, p_mpa: f64) -> f64 {
    mu_iapws2008(t_k, rho_if97(t_k, p_mpa)) * 1e3
}

// ----------------------------------------------------------------------------
// Bradley & Pitzer (1979) dielectric constant, Eqs. (1)-(4), Table I.
// T in K, P in bar.
// ----------------------------------------------------------------------------

const BP_U: [f64; 9] = [
    3.4279e2, -5.0866e-3, 9.4690e-7, -2.0525, 3.1159e3, -1.8289e2, -8.0325e3, 4.2142e6, 2.1417,
];

/// Static dielectric constant of water. Port of dielectric_constant.
pub fn dielectric_constant(t_k: f64, p_bar: f64) -> f64 {
    let u = &BP_U;
    let d1000 = u[0] * (u[1] * t_k + u[2] * t_k.powi(2)).exp();
    let c = u[3] + u[4] / (u[5] + t_k);
    let b = u[6] + u[7] / t_k + u[8] * t_k;
    d1000 + c * ((b + p_bar) / (b + 1000.0)).ln()
}
