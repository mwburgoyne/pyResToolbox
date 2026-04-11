use pyo3::prelude::*;
use crate::critical_properties;

// Constants matching pyrestoolbox/constants/constants.py
const R: f64 = 10.731577089016;
const MW_AIR: f64 = 28.97;
const DEGF2R: f64 = 459.67;

// BNS EOS component arrays (same as zfactor module)
const BNS_MWS: [f64; 5] = [44.01, 34.082, 28.014, 2.016, 0.0];
const BNS_TCS: [f64; 5] = [547.416, 672.120, 227.160, 47.430, 1.0];
const BNS_PCS: [f64; 5] = [1069.51, 1299.97, 492.84, 187.5300, 1.0];
const BNS_VCVIS: [f64; 5] = [1.46352, 1.46808, 1.35526, 0.68473, 0.0];

// LBC polynomial coefficients
const A_LBC: [f64; 5] = [0.1023, 0.023364, 0.058533, -0.0392852, 0.00926279];

// =========================================================================
// Lee-Gonzalez-Eakin viscosity (conventional method)
// =========================================================================

/// LGE viscosity for a single pressure point.
/// p_psia: pressure, t_degR: temperature in Rankine, sg: gas SG, zee: Z-factor
/// Returns viscosity in cP.
pub fn lge_viscosity(p_psia: f64, t_degr: f64, sg: f64, zee: f64) -> f64 {
    let m = MW_AIR * sg;
    let rho = m * p_psia / (t_degr * zee * R * 62.37);
    let b = 3.448 + (986.4 / t_degr) + (0.01009 * m);
    let c = 2.447 - (0.2224 * b);
    let a = (9.379 + (0.01607 * m)) * t_degr.powf(1.5) / (209.2 + (19.26 * m) + t_degr);
    a * 0.0001 * (b * rho.powf(c)).exp()
}

/// LGE viscosity exposed to Python (scalar).
#[pyfunction]
pub fn gas_ug_lge(
    p_psia: f64,
    sg: f64,
    degf: f64,
    zee: f64,
) -> PyResult<f64> {
    let t_degr = degf + DEGF2R;
    Ok(lge_viscosity(p_psia, t_degr, sg, zee))
}

// =========================================================================
// LBC viscosity (BNS method)
// =========================================================================

/// Precomputed LBC mixture parameters that are constant for a given composition and temperature.
pub struct LbcParams {
    pub u0: f64,        // Dilute gas mixture viscosity (Herning-Zippener)
    pub eta_mix: f64,   // LBC reducing parameter
    pub rhoc: f64,      // Critical density (1/Vc_mix)
}

/// Compute LBC parameters that don't depend on pressure.
/// These can be reused across multiple pressure evaluations.
pub fn lbc_params(
    degf: f64,
    sg: f64,
    co2: f64, h2s: f64, n2: f64, h2: f64,
) -> LbcParams {
    let deg_r = degf + DEGF2R;
    let zi = [co2, h2s, n2, h2, 1.0 - co2 - h2s - n2 - h2];

    let mut mws = BNS_MWS;
    let mut tcs = BNS_TCS;
    let mut pcs = BNS_PCS;
    let mut vcvis = BNS_VCVIS;

    // Compute HC properties
    let inert_sum = co2 + h2s + n2 + h2;
    let sg_hc = if inert_sum < 1.0 {
        let raw = (sg - (co2 * mws[0] + h2s * mws[1] + n2 * mws[2] + h2 * mws[3]) / MW_AIR) / (1.0 - inert_sum);
        if raw < 0.553779772 { 0.553779772 } else { raw }
    } else {
        0.75
    };
    let hc_gas_mw = sg_hc * MW_AIR;

    mws[4] = hc_gas_mw;
    let (tpc_hc, ppc_hc, _) = critical_properties::bns_pseudocritical_internal(
        hc_gas_mw / MW_AIR, 0.0, 0.0, 0.0, 0.0
    );
    tcs[4] = tpc_hc;
    pcs[4] = ppc_hc;
    vcvis[4] = 0.0576710 * (hc_gas_mw - 16.0425) + 1.44383;

    // Stiel-Thodos dilute gas viscosity per component
    let mut ui = [0.0; 5];
    for i in 0..5 {
        let tr = deg_r / tcs[i];
        let tc_k = tcs[i] * 5.0 / 9.0;
        let pc_atm = pcs[i] / 14.696;
        let eta_st = tc_k.powf(1.0 / 6.0) / (mws[i].sqrt() * pc_atm.powf(2.0 / 3.0));
        if tr <= 1.5 {
            ui[i] = 34e-5 * tr.powf(0.94) / eta_st;
        } else {
            let arg = (4.58 * tr - 1.67).max(1e-30);
            ui[i] = 17.78e-5 * arg.powf(5.0 / 8.0) / eta_st;
        }
    }

    // Herning-Zippener
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..5 {
        let sqrt_mw = mws[i].sqrt();
        num += zi[i] * ui[i] * sqrt_mw;
        den += zi[i] * sqrt_mw;
    }
    let u0 = num / den;

    // LBC mixture parameters
    let mut vc_sum = 0.0;
    let mut tc_k_sum = 0.0;
    let mut mw_sum = 0.0;
    let mut pc_atm_sum = 0.0;
    for i in 0..5 {
        vc_sum += vcvis[i] * zi[i];
        tc_k_sum += zi[i] * tcs[i] * 5.0 / 9.0;
        mw_sum += zi[i] * mws[i];
        pc_atm_sum += zi[i] * pcs[i] / 14.696;
    }
    let rhoc = 1.0 / vc_sum;
    let eta_mix = tc_k_sum.abs().powf(1.0 / 6.0)
        / (mw_sum.abs().sqrt() * pc_atm_sum.abs().powf(2.0 / 3.0));

    LbcParams { u0, eta_mix, rhoc }
}

/// LBC viscosity at a single pressure point given precomputed params.
pub fn lbc_viscosity_with_params(
    p_psia: f64,
    deg_r: f64,
    zee: f64,
    params: &LbcParams,
) -> f64 {
    let rhor = p_psia / (zee * R * deg_r * params.rhoc);
    let lhs = A_LBC[0] + rhor * (A_LBC[1] + rhor * (A_LBC[2] + rhor * (A_LBC[3] + rhor * A_LBC[4])));
    let lhs4 = lhs * lhs * lhs * lhs;
    (lhs4 - 1e-4) / params.eta_mix + params.u0
}

/// LBC viscosity exposed to Python (scalar, full computation).
#[pyfunction]
pub fn gas_ug_lbc(
    p_psia: f64,
    sg: f64,
    degf: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
    zee: f64,
) -> PyResult<f64> {
    let deg_r = degf + DEGF2R;
    let params = lbc_params(degf, sg, co2, h2s, n2, h2);
    Ok(lbc_viscosity_with_params(p_psia, deg_r, zee, &params))
}

// =========================================================================
// Batch (vectorized) viscosity functions
// =========================================================================

/// LGE viscosity for a batch of (pressure, z-factor) pairs.
#[pyfunction]
pub fn gas_ug_lge_batch(
    pressures: Vec<f64>,
    z_factors: Vec<f64>,
    sg: f64,
    degf: f64,
) -> PyResult<Vec<f64>> {
    let t_degr = degf + DEGF2R;
    let m = MW_AIR * sg;

    // Precompute temperature-dependent LGE coefficients
    let b = 3.448 + (986.4 / t_degr) + (0.01009 * m);
    let c = 2.447 - (0.2224 * b);
    let a = (9.379 + (0.01607 * m)) * t_degr.powf(1.5) / (209.2 + (19.26 * m) + t_degr);

    let result: Vec<f64> = pressures.iter().zip(z_factors.iter()).map(|(&p, &z)| {
        let rho = m * p / (t_degr * z * R * 62.37);
        a * 0.0001 * (b * rho.powf(c)).exp()
    }).collect();

    Ok(result)
}

/// LBC viscosity for a batch of (pressure, z-factor) pairs.
/// Precomputes LBC mixture parameters (u0, eta_mix, rhoc) once.
#[pyfunction]
pub fn gas_ug_lbc_batch(
    pressures: Vec<f64>,
    z_factors: Vec<f64>,
    sg: f64,
    degf: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
) -> PyResult<Vec<f64>> {
    let deg_r = degf + DEGF2R;
    let params = lbc_params(degf, sg, co2, h2s, n2, h2);

    let result: Vec<f64> = pressures.iter().zip(z_factors.iter()).map(|(&p, &z)| {
        lbc_viscosity_with_params(p, deg_r, z, &params)
    }).collect();

    Ok(result)
}
