use pyo3::prelude::*;

mod zfactor;
mod critical_properties;
mod gas_viscosity;
mod pseudopressure;
mod vlp;
mod dca;
mod matbal;
mod oil;
mod spycher_pruess;
mod vle;
mod bessel;
mod gwr;

/// Smoke test function called during import-time probe.
#[pyfunction]
fn _smoke_test() -> PyResult<()> {
    let _x: f64 = (2.0_f64).sqrt();
    Ok(())
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_smoke_test, m)?)?;

    // Critical properties
    m.add_function(wrap_pyfunction!(critical_properties::sutton_pseudocritical, m)?)?;
    m.add_function(wrap_pyfunction!(critical_properties::wichert_aziz_correction, m)?)?;
    m.add_function(wrap_pyfunction!(critical_properties::bns_pseudocritical, m)?)?;

    // Z-factor correlations (reduced-property input)
    m.add_function(wrap_pyfunction!(zfactor::dak_zfactor, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::hall_yarborough_zfactor, m)?)?;

    // Full pipeline functions
    m.add_function(wrap_pyfunction!(zfactor::dak_zfactor_full, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::hall_yarborough_zfactor_full, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::bns_zfactor_full, m)?)?;

    // Z-factor batch functions
    m.add_function(wrap_pyfunction!(zfactor::dak_zfactor_batch, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::hy_zfactor_batch, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::bns_zfactor_batch, m)?)?;

    // Gas viscosity
    m.add_function(wrap_pyfunction!(gas_viscosity::gas_ug_lge, m)?)?;
    m.add_function(wrap_pyfunction!(gas_viscosity::gas_ug_lbc, m)?)?;

    // Gas viscosity batch functions
    m.add_function(wrap_pyfunction!(gas_viscosity::gas_ug_lge_batch, m)?)?;
    m.add_function(wrap_pyfunction!(gas_viscosity::gas_ug_lbc_batch, m)?)?;

    // Pseudopressure integration
    m.add_function(wrap_pyfunction!(pseudopressure::gas_dmp_rust, m)?)?;

    // P/Z -> P solver
    m.add_function(wrap_pyfunction!(pseudopressure::gas_ponz2p_rust, m)?)?;

    // VLP segment loops
    m.add_function(wrap_pyfunction!(vlp::hb_fbhp_gas_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::hb_fbhp_oil_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::wg_fbhp_gas_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::wg_fbhp_oil_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::gray_fbhp_gas_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::gray_fbhp_oil_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::bb_fbhp_gas_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::bb_fbhp_oil_rust, m)?)?;

    // DCA hyperbolic grid search
    m.add_function(wrap_pyfunction!(dca::fit_hyperbolic_rust, m)?)?;
    m.add_function(wrap_pyfunction!(dca::fit_hyperbolic_cum_rust, m)?)?;

    // Material balance objective
    m.add_function(wrap_pyfunction!(matbal::matbal_objective_rust, m)?)?;

    // Oil density chain
    m.add_function(wrap_pyfunction!(oil::oil_deno_mccain_rust, m)?)?;
    m.add_function(wrap_pyfunction!(oil::oil_bo_mccain_rust, m)?)?;

    // Spycher-Pruess CO2-Brine solubility
    m.add_function(wrap_pyfunction!(spycher_pruess::co2_brine_solubility_rust, m)?)?;

    // VLE Flash Engine
    m.add_function(wrap_pyfunction!(vle::flash_tp_rust, m)?)?;
    m.add_function(wrap_pyfunction!(vle::calc_equilibrium_rust, m)?)?;

    // GWR inverse Laplace transform / influence tables
    m.add_function(wrap_pyfunction!(gwr::influence_tables_rust, m)?)?;
    m.add_function(wrap_pyfunction!(gwr::gwr_rust, m)?)?;

    Ok(())
}
