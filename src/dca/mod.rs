pub mod ransac;
pub mod hyperbolic;

use pyo3::prelude::*;
use hyperbolic::{fit_hyperbolic, fit_hyperbolic_cum};

/// Fit hyperbolic decline to rate-vs-time data using grid search + RANSAC.
///
/// Returns (qi, di, b, R2).
#[pyfunction]
pub fn fit_hyperbolic_rust(t: Vec<f64>, q: Vec<f64>) -> PyResult<(f64, f64, f64, f64)> {
    if t.len() != q.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "t and q must have the same length",
        ));
    }
    if t.len() < 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Need at least 3 data points",
        ));
    }
    match fit_hyperbolic(&t, &q) {
        Some(result) => Ok(result),
        None => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Hyperbolic fit failed: no valid b produced a positive qi and di",
        )),
    }
}

/// Fit hyperbolic decline to cumulative-production-vs-rate data using grid search + RANSAC.
///
/// Returns (qi, di, b, R2).
#[pyfunction]
pub fn fit_hyperbolic_cum_rust(np_cum: Vec<f64>, q: Vec<f64>) -> PyResult<(f64, f64, f64, f64)> {
    if np_cum.len() != q.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "np_cum and q must have the same length",
        ));
    }
    if np_cum.len() < 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Need at least 3 data points",
        ));
    }
    match fit_hyperbolic_cum(&np_cum, &q) {
        Some(result) => Ok(result),
        None => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Hyperbolic cumulative fit failed: no valid b produced a positive qi and di",
        )),
    }
}
