use pyo3::prelude::*;
use crate::critical_properties;
use crate::zfactor;
use crate::gas_viscosity;

const DEGF2R: f64 = 459.67;
#[allow(dead_code)] const MW_AIR: f64 = 28.97;
#[allow(dead_code)] const R: f64 = 10.731577089016;

// Gauss-Legendre 7-point nodes and weights
const GL7_NODES: [f64; 7] = [
    -0.9491079123427585,
    -0.7415311855993945,
    -0.4058451513773972,
    0.0,
    0.4058451513773972,
    0.7415311855993945,
    0.9491079123427585,
];
const GL7_WEIGHTS: [f64; 7] = [
    0.1294849661688697,
    0.2797053914892766,
    0.3818300505051189,
    0.4179591836734694,
    0.3818300505051189,
    0.2797053914892766,
    0.1294849661688697,
];

// Gauss-Legendre 10-point nodes and weights
const GL10_NODES: [f64; 10] = [
    -0.9739065285171717,
    -0.8650633666889845,
    -0.6794095682990244,
    -0.4333953941292472,
    -0.1488743389816312,
    0.1488743389816312,
    0.4333953941292472,
    0.6794095682990244,
    0.8650633666889845,
    0.9739065285171717,
];
const GL10_WEIGHTS: [f64; 10] = [
    0.0666713443086881,
    0.1494513491505806,
    0.2190863625159820,
    0.2692667193099963,
    0.2955242247147529,
    0.2955242247147529,
    0.2692667193099963,
    0.2190863625159820,
    0.1494513491505806,
    0.0666713443086881,
];

/// Method selector for the full pipeline
enum ZMethod {
    DakSut,
    HySut,
    Bns,
}

/// Evaluate Z-factor at a single pressure using the selected method.
/// All critical property computation is internal — no Python round-trips.
fn eval_z(
    p_psia: f64,
    method: &ZMethod,
    deg_r: f64,
    _degf: f64,
    // Sutton path params (precomputed)
    tpc_sut: f64,
    ppc_sut: f64,
    // BNS path params (precomputed)
    tpc_bns: f64,
    ppc_bns: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
) -> f64 {
    match method {
        ZMethod::DakSut => {
            let pr = p_psia / ppc_sut;
            let tr = deg_r / tpc_sut;
            zfactor::dak_core_pub(pr, tr)
        }
        ZMethod::HySut => {
            let pr = p_psia / ppc_sut;
            let tr = deg_r / tpc_sut;
            zfactor::hy_core_pub(pr, tr)
        }
        ZMethod::Bns => {
            zfactor::bns_zfactor_core_pub(
                p_psia, deg_r, co2, h2s, n2, h2, tpc_bns, ppc_bns
            )
        }
    }
}

/// Evaluate viscosity at a single pressure using the appropriate method.
fn eval_ug(
    p_psia: f64,
    deg_r: f64,
    sg: f64,
    zee: f64,
    method: &ZMethod,
    lbc_params: &Option<gas_viscosity::LbcParams>,
) -> f64 {
    match method {
        ZMethod::Bns => {
            if let Some(ref params) = lbc_params {
                gas_viscosity::lbc_viscosity_with_params(p_psia, deg_r, zee, params)
            } else {
                gas_viscosity::lge_viscosity(p_psia, deg_r, sg, zee)
            }
        }
        _ => {
            gas_viscosity::lge_viscosity(p_psia, deg_r, sg, zee)
        }
    }
}

/// Gauss-Legendre integration of 2p/(mu*Z) over [lo, hi] using n-point rule.
/// Everything stays in Rust — Z-factor, viscosity, quadrature.
fn gl_integrate(
    lo: f64,
    hi: f64,
    nodes: &[f64],
    weights: &[f64],
    method: &ZMethod,
    deg_r: f64,
    degf: f64,
    sg: f64,
    tpc_sut: f64,
    ppc_sut: f64,
    tpc_bns: f64,
    ppc_bns: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
    lbc_params: &Option<gas_viscosity::LbcParams>,
) -> f64 {
    let p_mid = (lo + hi) * 0.5;
    let p_half = (hi - lo) * 0.5;
    let mut result = 0.0;

    for (i, &node) in nodes.iter().enumerate() {
        let p_eval = p_mid + p_half * node;
        let zee = eval_z(p_eval, method, deg_r, degf, tpc_sut, ppc_sut, tpc_bns, ppc_bns, co2, h2s, n2, h2);
        let ug = eval_ug(p_eval, deg_r, sg, zee, method, lbc_params);
        let mugz = ug * zee;
        result += weights[i] * 2.0 * p_eval / mugz;
    }

    p_half * result
}

/// Full pseudopressure integration: gas_dmp entirely in Rust.
/// zmethod: "DAK", "HY", or "BNS"
/// cmethod: "SUT" or "BNS" (PMC not accelerated — falls back to Python)
#[pyfunction]
pub fn gas_dmp_rust(
    p1: f64,
    p2: f64,
    degf: f64,
    sg: f64,
    zmethod: &str,
    cmethod: &str,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
    tc: f64,
    pc: f64,
) -> PyResult<f64> {
    if p1 == p2 {
        return Ok(0.0);
    }

    let deg_r = degf + DEGF2R;

    // Determine method
    let method = match (zmethod, cmethod) {
        ("DAK", "SUT") | ("DAK", _) if cmethod == "SUT" => ZMethod::DakSut,
        ("HY", "SUT") | ("HY", _) if cmethod == "SUT" => ZMethod::HySut,
        ("BNS", _) | ("BUR", _) => ZMethod::Bns,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                format!("Rust acceleration not available for zmethod={}, cmethod={}", zmethod, cmethod)
            ));
        }
    };

    // Precompute critical properties
    let (tpc_sut, ppc_sut) = if tc > 0.0 && pc > 0.0 {
        (tc, pc)
    } else {
        match method {
            ZMethod::DakSut | ZMethod::HySut => {
                critical_properties::sutton_wa_internal(sg, co2, h2s, n2)
                    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?
            }
            ZMethod::Bns => (0.0, 0.0), // Not used
        }
    };

    let (tpc_bns, ppc_bns) = match method {
        ZMethod::Bns => {
            let (t, p, _) = critical_properties::bns_pseudocritical_internal(sg, co2, h2s, n2, h2);
            (t, p)
        }
        _ => (0.0, 0.0),
    };

    // Precompute LBC params for BNS method
    let lbc_p = match method {
        ZMethod::Bns => Some(gas_viscosity::lbc_params(degf, sg, co2, h2s, n2, h2)),
        _ => None,
    };

    // Two-tier Gauss-Legendre integration
    let result_7 = gl_integrate(
        p1, p2, &GL7_NODES, &GL7_WEIGHTS,
        &method, deg_r, degf, sg,
        tpc_sut, ppc_sut, tpc_bns, ppc_bns,
        co2, h2s, n2, h2, &lbc_p,
    );
    let result_10 = gl_integrate(
        p1, p2, &GL10_NODES, &GL10_WEIGHTS,
        &method, deg_r, degf, sg,
        tpc_sut, ppc_sut, tpc_bns, ppc_bns,
        co2, h2s, n2, h2, &lbc_p,
    );

    if result_10.abs() < 1e-30 || (result_10 - result_7).abs() / result_10.abs() < 1e-5 {
        Ok(result_10)
    } else {
        // Split into two subintervals, integrate each with 10-point
        let p_mid = (p1 + p2) * 0.5;
        let r_lo = gl_integrate(
            p1, p_mid, &GL10_NODES, &GL10_WEIGHTS,
            &method, deg_r, degf, sg,
            tpc_sut, ppc_sut, tpc_bns, ppc_bns,
            co2, h2s, n2, h2, &lbc_p,
        );
        let r_hi = gl_integrate(
            p_mid, p2, &GL10_NODES, &GL10_WEIGHTS,
            &method, deg_r, degf, sg,
            tpc_sut, ppc_sut, tpc_bns, ppc_bns,
            co2, h2s, n2, h2, &lbc_p,
        );
        Ok(r_lo + r_hi)
    }
}
