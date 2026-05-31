/// VLE Flash Engine for Multi-Gas Phase Equilibria (Rust acceleration).
///
/// Implements the Soreide-Whitson Peng-Robinson EOS-based flash calculation
/// with up to 16 components for gas-saturated brine.
///
/// Exports:
///   flash_tp_rust — Single TP flash (AQ or NA mode)

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

/// Single TP flash for the S&W VLE engine.
///
/// # Arguments
/// * `t_k` - Temperature in Kelvin
/// * `p_pa` - Pressure in Pascal
/// * `z` - Feed composition (will be normalized)
/// * `comp_names` - Component names (e.g., ["H2O", "CH4", "CO2"])
/// * `salinity` - NaCl molality (mol/kg water). Informational only; the caller
///                is responsible for baking salinity into `gamma`.
/// * `mode` - "AQ" for aqueous phase BIPs, "NA" for non-aqueous
/// * `gamma` - Activity coefficient array (same length as z). Must be computed
///             by the caller. For framework='proposed', specialized ks models
///             (Dubessy/Akinfiev/Li/Mao-Duan/Duan-Sun) live on the Python side
///             — this Rust entry point must not substitute its own S&W Eq 8
///             fallback silently. Use `1.0` for gases in freshwater.
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
    mode: String,
    gamma: Vec<f64>,
) -> PyResult<(f64, Vec<f64>, Vec<f64>, bool)> {
    // Salinity arg is retained for signature stability / future diagnostics;
    // currently unused because gamma is always caller-supplied.
    let _ = salinity;

    // Validate inputs
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

    let mode_aq = match mode.as_str() {
        "AQ" => true,
        "NA" => false,
        _ => {
            return Err(PyValueError::new_err(
                "mode must be 'AQ' or 'NA'",
            ));
        }
    };

    let comp_indices = resolve_comp_names(&comp_names)?;

    let (v, x, y, converged) = flash::flash_tp(
        t_k,
        p_pa,
        &z,
        &comp_indices,
        mode_aq,
        &gamma,
        200,
        1e-10,
    );

    Ok((v, x, y, converged))
}

