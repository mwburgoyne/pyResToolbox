//! Gas-free NaCl brine viscosity: IAPWS-2008 water (IF97 density) times the
//! ion-additive Jones-Dole salt ratio times Kestin's pressure factor.
//!
//! Direct port of brine/viscosity_route.py brine_viscosity(T_K, P_MPa, m=m)
//! on its default route ('jones_dole', pressure_term=True), restricted to a
//! pure NaCl brine. Used by the VLP segment marches for the water phase.
//! The Python chain is the reference; agreement is checked to 1e-9 relative
//! in test_rust_acceleration.py.
//!
//! No range checks: callers keep T within 273.15-623.15 K and P within
//! 0-100 MPa (IF97 Region 1), as the Python chain enforces by raising.

use pyo3::prelude::*;

pub mod salt;
pub mod water;

/// Gas-free NaCl brine viscosity in mPa s (cP). T in K, P in MPa, m in mol/kg
/// water. m <= 0 returns pure-water viscosity.
pub fn brine_viscosity_nacl(t_k: f64, p_mpa: f64, m: f64) -> f64 {
    let mu_w = water::mu_water(t_k, p_mpa);
    if m <= 0.0 {
        return mu_w;
    }
    let ratio = salt::jones_dole_ratio(t_k, p_mpa, m, mu_w)
        * salt::kestin_pressure_factor(t_k, p_mpa, m);
    mu_w * ratio
}

/// Python-facing wrapper for parity testing against
/// viscosity_route.brine_viscosity(T_K, P_MPa, m=m).
#[pyfunction]
pub fn brine_viscosity_nacl_rust(t_k: f64, p_mpa: f64, m: f64) -> f64 {
    brine_viscosity_nacl(t_k, p_mpa, m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_water_matches_iapws() {
        // Huber et al. (2009) Table 4 check value: 298.15 K, 998 kg/m3.
        let mu = water::mu_iapws2008(298.15, 998.0);
        assert!((mu * 1e6 - 889.735100).abs() < 1e-5, "{mu}");
    }

    #[test]
    fn if97_density_check_value() {
        // IF97 Table 5: 300 K, 3 MPa -> v = 0.100215168e-2 m3/kg.
        let rho = water::rho_if97(300.0, 3.0);
        assert!((1.0 / rho - 0.100215168e-2).abs() < 1e-11, "{rho}");
    }

    #[test]
    fn salt_raises_viscosity() {
        let fresh = brine_viscosity_nacl(350.0, 20.0, 0.0);
        let salty = brine_viscosity_nacl(350.0, 20.0, 2.0);
        assert!(salty > fresh * 1.1 && salty < fresh * 1.5);
    }
}
