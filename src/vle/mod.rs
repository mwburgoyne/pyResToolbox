//! VLE Flash Engine for Multi-Gas Phase Equilibria (Rust acceleration).
//!
//! Implements the Soreide-Whitson Peng-Robinson EOS flash on S&W's own
//! "procedure 1" (S&W 1992 Eqs 10a/10b: kij_AQ in the aqueous phase, kij_NA
//! in the non-aqueous phase) for gas-saturated brine, 'default' framework.
//!
//! Exports:
//!   flash_tp_rust - single procedure-1 TP flash

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

pub mod components;
pub mod alpha;
pub mod bip;
pub mod eos;
pub mod fugacity;
pub mod rachford_rice;
pub mod k_init;
pub mod flash;

use components::comp_index;

/// Resolve component names to indices. Returns error if any name is unknown
/// or H2O is missing.
fn resolve_comp_names(comp_names: &[String]) -> PyResult<Vec<usize>> {
    let mut indices = Vec::with_capacity(comp_names.len());
    let mut has_water = false;
    for name in comp_names {
        match comp_index(name.as_str()) {
            Some(idx) => {
                if idx == components::IDX_H2O {
                    has_water = true;
                }
                indices.push(idx);
            }
            None => {
                return Err(PyValueError::new_err(format!(
                    "Unknown component: {}. Supported: {:?}",
                    name,
                    components::COMPONENT_NAMES
                )));
            }
        }
    }
    if !has_water {
        return Err(PyValueError::new_err("H2O must be in component list"));
    }
    Ok(indices)
}

/// Single procedure-1 TP flash for the S&W VLE engine ('default' framework:
/// S&W water alpha, freshwater kij_AQ plus embedded salinity delta_kij in the
/// aqueous phase, kij_NA in the non-aqueous phase).
///
/// # Arguments
/// * `t_k` - Temperature in Kelvin
/// * `p_pa` - Pressure in Pascal
/// * `z` - Feed composition (will be normalized)
/// * `comp_names` - Component names (e.g., ["H2O", "CH4", "CO2"])
/// * `salinity` - NaCl molality (mol/kg water). Enters the water alpha and the
///                embedded delta_kij on kij_AQ.
/// * `gamma` - Aqueous activity coefficients (same length as z), computed by
///             the caller; all-ones under embedded salinity. This entry point
///             never substitutes its own Sechenov values.
///
/// # Returns
/// Tuple (V, x_vec, y_vec, converged) where V is vapor fraction, x is liquid
/// comp, y is vapor comp, and converged flags whether the SS loop converged.
#[pyfunction]
pub fn flash_tp_rust(
    t_k: f64,
    p_pa: f64,
    z: Vec<f64>,
    comp_names: Vec<String>,
    salinity: f64,
    gamma: Vec<f64>,
) -> PyResult<(f64, Vec<f64>, Vec<f64>, bool)> {
    if !t_k.is_finite() || t_k <= 0.0 {
        return Err(PyValueError::new_err("Temperature must be positive and finite"));
    }
    if !p_pa.is_finite() || p_pa <= 0.0 {
        return Err(PyValueError::new_err("Pressure must be positive and finite"));
    }
    if z.len() != comp_names.len() {
        return Err(PyValueError::new_err(
            "z and comp_names must have the same length",
        ));
    }
    if gamma.len() != comp_names.len() {
        return Err(PyValueError::new_err(
            "gamma and comp_names must have the same length",
        ));
    }

    let comp_indices = resolve_comp_names(&comp_names)?;

    Ok(flash::flash_tp(t_k, p_pa, &z, &comp_indices, &gamma, salinity, 200, 1e-10))
}
