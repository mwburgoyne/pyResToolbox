use pyo3::prelude::*;

mod solubility;

/// Spycher-Pruess CO2-Brine mutual solubility calculation.
///
/// Computes the phase-partitioning of CO2 and H2O between an aqueous brine
/// phase and a CO2-rich gas phase at elevated temperatures and pressures.
///
/// # Arguments
/// * `p_bar`  - Pressure in bar
/// * `deg_c`  - Temperature in degrees Celsius
/// * `ppm`    - NaCl concentration in weight parts per million of brine
///
/// # Returns
/// Tuple of (x_co2, y_co2, y_h2o, rho_gas, gas_z) where:
/// * `x_co2`   - Mole fraction of CO2 in the aqueous phase
/// * `y_co2`   - Mole fraction of CO2 in the gas phase
/// * `y_h2o`   - Mole fraction of H2O in the gas phase
/// * `rho_gas`  - CO2-rich gas density (g/cm3)
/// * `gas_z`   - Gas phase Z-factor
#[pyfunction]
pub fn co2_brine_solubility_rust(
    p_bar: f64,
    deg_c: f64,
    ppm: f64,
) -> PyResult<(f64, f64, f64, f64, f64)> {
    Ok(solubility::co2_brine_solubility(p_bar, deg_c, ppm))
}
