use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

// Constants must match pyrestoolbox/constants/constants.py exactly
const R: f64 = 10.731577089016;
const MW_AIR: f64 = 28.97;
const MW_CO2: f64 = 44.01;
const MW_H2S: f64 = 34.082;
const MW_N2: f64 = 28.014;
const MW_H2: f64 = 2.016;

/// Sutton (1985) pseudocritical properties from gas specific gravity.
/// Returns (Ppc_psia, Tpc_degR) for the HC+impurity mixture using
/// Kay's mixing rule for N2, CO2, H2S with Wichert-Aziz correction.
///
/// Note: Uses the same MW values as the Python implementation in gas.py
/// for the SG calculation (N2=28.01, CO2=44.01, H2S=34.1), which differ
/// slightly from the constants module values used in BNS.
#[pyfunction]
pub fn sutton_pseudocritical(
    sg: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
) -> PyResult<(f64, f64)> {
    let hc_frac = 1.0 - n2 - co2 - h2s;
    if hc_frac <= 0.0 {
        return Err(PyValueError::new_err(
            "SUT method requires hydrocarbon fraction > 0 (n2 + co2 + h2s must be < 1.0)"
        ));
    }
    // MW values match gas.py line 391: N2=28.01, CO2=44.01, H2S=34.1
    let sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / MW_AIR) / hc_frac;

    let ppc_hc = 756.8 - 131.0 * sg_hc - 3.6 * sg_hc * sg_hc;
    let tpc_hc = 169.2 + 349.5 * sg_hc - 74.0 * sg_hc * sg_hc;

    // Kay's mixing rule
    let ppc_star = hc_frac * ppc_hc + n2 * 507.5 + co2 * 1071.0 + h2s * 1306.0;
    let tpc_star = hc_frac * tpc_hc + n2 * 239.26 + co2 * 547.58 + h2s * 672.35;

    // Wichert-Aziz acid gas correction
    let a = co2 + h2s;
    let eps = 120.0 * (a.powf(0.9) - a.powf(1.6))
            + 15.0 * (h2s.powf(0.5) - h2s.powf(4.0));

    let tpc = tpc_star - eps;
    let ppc = ppc_star * (tpc_star - eps) / (tpc_star + h2s * (1.0 - h2s) * eps);

    Ok((ppc, tpc))
}

/// Wichert-Aziz (1972) acid gas correction applied to pseudocritical properties.
/// Takes uncorrected Ppc, Tpc and returns corrected values.
#[pyfunction]
pub fn wichert_aziz_correction(
    ppc: f64,
    tpc: f64,
    co2_frac: f64,
    h2s_frac: f64,
) -> PyResult<(f64, f64)> {
    let a = co2_frac + h2s_frac;
    let eps = 120.0 * (a.powf(0.9) - a.powf(1.6))
            + 15.0 * (h2s_frac.powf(0.5) - h2s_frac.powf(4.0));

    let tpc_corr = tpc - eps;
    let ppc_corr = ppc * (tpc - eps) / (tpc + h2s_frac * (1.0 - h2s_frac) * eps);

    Ok((ppc_corr, tpc_corr))
}

/// BNS pseudocritical properties using Burgoyne (2025) correlation.
/// Returns (Ppc_psia, Tpc_degR) for the hydrocarbon pseudo-component.
/// Does NOT apply Wichert-Aziz — BNS mixing rules handle acid gas internally.
#[pyfunction]
pub fn bns_pseudocritical(
    sg: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
) -> PyResult<(f64, f64)> {
    let inert_sum = co2 + h2s + n2 + h2;
    let sg_hc = if inert_sum < 1.0 {
        let raw = (sg - (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2) / MW_AIR) / (1.0 - inert_sum);
        if raw < 0.553779772 { 0.553779772 } else { raw }
    } else {
        0.75
    };

    let x = (MW_AIR * sg_hc - 16.0425_f64).max(0.0);

    // Gas Condensate regime (default)
    let a_gc: f64 = 1098.10948;
    let b_gc: f64 = 101.529237;
    let c_tc: f64 = 343.008;
    let tpc = a_gc * x / (b_gc + x) + c_tc;

    let vc_slope: f64 = 0.170931432;
    let vc_on_zc = vc_slope * x + 5.518525872412144;
    let ppc = R * tpc / vc_on_zc;

    Ok((ppc, tpc))
}

/// Internal: compute BNS pseudocritical properties (not exposed to Python).
/// Returns (tpc, ppc, sg_hc) for use by z_bur full pipeline.
pub fn bns_pseudocritical_internal(
    sg: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
) -> (f64, f64, f64) {
    let inert_sum = co2 + h2s + n2 + h2;
    let sg_hc = if inert_sum < 1.0 {
        let raw = (sg - (co2 * MW_CO2 + h2s * MW_H2S + n2 * MW_N2 + h2 * MW_H2) / MW_AIR) / (1.0 - inert_sum);
        if raw < 0.553779772 { 0.553779772 } else { raw }
    } else {
        0.75
    };

    let x = (MW_AIR * sg_hc - 16.0425_f64).max(0.0);

    let a_gc: f64 = 1098.10948;
    let b_gc: f64 = 101.529237;
    let c_tc: f64 = 343.008;
    let tpc = a_gc * x / (b_gc + x) + c_tc;

    let vc_slope: f64 = 0.170931432;
    let vc_on_zc = vc_slope * x + 5.518525872412144;
    let ppc = R * tpc / vc_on_zc;

    (tpc, ppc, sg_hc)
}

/// Internal: Sutton + Wichert-Aziz for DAK/HY full pipeline.
/// Returns (tpc_corrected, ppc_corrected).
pub fn sutton_wa_internal(
    sg: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
) -> Result<(f64, f64), String> {
    let hc_frac = 1.0 - n2 - co2 - h2s;
    if hc_frac <= 0.0 {
        return Err("SUT method requires hydrocarbon fraction > 0".to_string());
    }
    // MW values match gas.py line 391: N2=28.01, CO2=44.01, H2S=34.1
    let sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / MW_AIR) / hc_frac;

    let ppc_hc = 756.8 - 131.0 * sg_hc - 3.6 * sg_hc * sg_hc;
    let tpc_hc = 169.2 + 349.5 * sg_hc - 74.0 * sg_hc * sg_hc;

    let ppc_star = hc_frac * ppc_hc + n2 * 507.5 + co2 * 1071.0 + h2s * 1306.0;
    let tpc_star = hc_frac * tpc_hc + n2 * 239.26 + co2 * 547.58 + h2s * 672.35;

    let a = co2 + h2s;
    let eps = 120.0 * (a.powf(0.9) - a.powf(1.6))
            + 15.0 * (h2s.powf(0.5) - h2s.powf(4.0));

    let tpc = tpc_star - eps;
    let ppc = ppc_star * (tpc_star - eps) / (tpc_star + h2s * (1.0 - h2s) * eps);

    Ok((tpc, ppc))
}
