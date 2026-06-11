/// VLP segment loop acceleration module.
/// Exposes 8 pyfunction exports for the 4 VLP methods x 2 fluid types.
/// All entry points route through the shared segment march in march.rs.
/// A diverged march (pressure below 1 psia) maps to a Python RuntimeError
/// with the same message as the pure-Python implementation in nodal.py.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

pub mod constants;
pub mod pvt_helpers;
pub mod ift;
pub mod friction;
pub mod holdup_wg;
pub mod holdup_gray;
pub mod holdup_bb;
pub mod march;
pub mod static_column;

/// HB gas VLP segment loop.
#[pyfunction]
pub fn hb_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    march::hb_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    )
    .map_err(PyRuntimeError::new_err)
}

/// HB oil VLP segment loop.
#[pyfunction]
pub fn hb_fbhp_oil_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qt_stbpd: f64, gor: f64, wc: f64, pb: f64, rsb: f64, sgsp: f64,
    rsb_scale: f64, injection: bool, theta: f64,
    vis_frac: f64, rsb_frac: f64,
) -> PyResult<f64> {
    march::hb_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    )
    .map_err(PyRuntimeError::new_err)
}

/// WG gas VLP segment loop.
#[pyfunction]
pub fn wg_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    march::wg_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    )
    .map_err(PyRuntimeError::new_err)
}

/// WG oil VLP segment loop.
#[pyfunction]
pub fn wg_fbhp_oil_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qt_stbpd: f64, gor: f64, wc: f64, pb: f64, rsb: f64, sgsp: f64,
    rsb_scale: f64, injection: bool, theta: f64,
    vis_frac: f64, rsb_frac: f64,
) -> PyResult<f64> {
    march::wg_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    )
    .map_err(PyRuntimeError::new_err)
}

/// Gray gas VLP segment loop.
#[pyfunction]
pub fn gray_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    march::gray_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    )
    .map_err(PyRuntimeError::new_err)
}

/// Gray oil VLP segment loop.
#[pyfunction]
pub fn gray_fbhp_oil_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qt_stbpd: f64, gor: f64, wc: f64, pb: f64, rsb: f64, sgsp: f64,
    rsb_scale: f64, injection: bool, theta: f64,
    vis_frac: f64, rsb_frac: f64,
) -> PyResult<f64> {
    march::gray_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    )
    .map_err(PyRuntimeError::new_err)
}

/// BB gas VLP segment loop.
#[pyfunction]
pub fn bb_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    march::bb_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    )
    .map_err(PyRuntimeError::new_err)
}

/// BB oil VLP segment loop.
#[pyfunction]
pub fn bb_fbhp_oil_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qt_stbpd: f64, gor: f64, wc: f64, pb: f64, rsb: f64, sgsp: f64,
    rsb_scale: f64, injection: bool, theta: f64,
    vis_frac: f64, rsb_frac: f64,
) -> PyResult<f64> {
    march::bb_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    )
    .map_err(PyRuntimeError::new_err)
}
