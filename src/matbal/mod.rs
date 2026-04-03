use pyo3::prelude::*;

mod objective;

/// Material balance objective function for L-BFGS-B optimization.
/// Computes coefficient of variation of OOIP estimates.
#[pyfunction]
pub fn matbal_objective_rust(
    params: Vec<f64>,
    param_names: Vec<String>,
    f_arr: Vec<f64>,
    eo_arr: Vec<f64>,
    eg_arr: Vec<f64>,
    p_field: Vec<f64>,
    boi: f64,
    base_m: f64,
    base_cf: f64,
    base_cw: f64,
    base_sw_i: f64,
) -> PyResult<f64> {
    Ok(objective::matbal_objective(
        &params, &param_names, &f_arr, &eo_arr, &eg_arr, &p_field,
        boi, base_m, base_cf, base_cw, base_sw_i,
    ))
}
