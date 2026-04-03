/// VLP segment loop acceleration module.
/// Exposes 8 pyfunction exports for the 4 VLP methods × 2 fluid types.

use pyo3::prelude::*;

pub mod pvt_helpers;
pub mod ift;
pub mod friction;
pub mod holdup_hb;
pub mod holdup_wg;
pub mod holdup_gray;
pub mod holdup_bb;
pub mod segment_gas;
pub mod segment_oil;
pub mod static_column;

/// HB gas VLP segment loop.
#[pyfunction]
pub fn hb_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    Ok(segment_gas::hb_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    ))
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
    Ok(segment_oil::hb_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    ))
}

/// WG gas VLP segment loop.
#[pyfunction]
pub fn wg_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    Ok(segment_gas::wg_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    ))
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
    Ok(segment_oil::wg_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    ))
}

/// Gray gas VLP segment loop.
#[pyfunction]
pub fn gray_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    Ok(segment_gas::gray_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    ))
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
    Ok(segment_oil::gray_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    ))
}

/// BB gas VLP segment loop.
#[pyfunction]
pub fn bb_fbhp_gas_rust(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> PyResult<f64> {
    Ok(segment_gas::bb_fbhp_gas(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta,
    ))
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
    Ok(segment_oil::bb_fbhp_oil(
        thp, api, gsg, tid, rough, length, tht, bht, wsg,
        qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
        vis_frac, rsb_frac,
    ))
}
