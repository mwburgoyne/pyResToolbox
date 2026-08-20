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
    0.219086362515982,
    0.2692667193099963,
    0.2955242247147529,
    0.2955242247147529,
    0.2692667193099963,
    0.219086362515982,
    0.1494513491505806,
    0.0666713443086881,
];

/// Method selector for the full pipeline
enum ZMethod {
    DakSut,
    HySut,
    Bns,
}

/// Gas state at a fixed temperature and composition: everything the Z-factor
/// and viscosity evaluators need, resolved once at the entry point.
///
/// These eleven values used to be threaded through eval_z, gl_integrate,
/// p_over_z and solve_ponz2p_single one argument at a time, and repeated
/// verbatim at every call site.
struct GasState<'a> {
    method: ZMethod,
    deg_r: f64,
    sg: f64,
    // Sutton path (precomputed)
    tpc_sut: f64,
    ppc_sut: f64,
    // BNS path (precomputed)
    tpc_bns: f64,
    ppc_bns: f64,
    co2: f64,
    h2s: f64,
    n2: f64,
    h2: f64,
    lbc_params: &'a Option<gas_viscosity::LbcParams>,
}

impl GasState<'_> {
    /// Z-factor at a single pressure. All critical property computation is
    /// internal — no Python round-trips.
    fn eval_z(&self, p_psia: f64) -> f64 {
        match self.method {
            ZMethod::DakSut => {
                zfactor::dak_core_pub(p_psia / self.ppc_sut, self.deg_r / self.tpc_sut)
            }
            ZMethod::HySut => {
                zfactor::hy_core_pub(p_psia / self.ppc_sut, self.deg_r / self.tpc_sut)
            }
            ZMethod::Bns => zfactor::bns_zfactor_core_pub(
                p_psia, self.deg_r, self.co2, self.h2s, self.n2, self.h2,
                self.tpc_bns, self.ppc_bns,
            ),
        }
    }

    /// Viscosity at a single pressure, using LBC only on the BNS path.
    fn eval_ug(&self, p_psia: f64, zee: f64) -> f64 {
        match self.method {
            ZMethod::Bns => match self.lbc_params {
                Some(params) => {
                    gas_viscosity::lbc_viscosity_with_params(p_psia, self.deg_r, zee, params)
                }
                None => gas_viscosity::lge_viscosity(p_psia, self.deg_r, self.sg, zee),
            },
            _ => gas_viscosity::lge_viscosity(p_psia, self.deg_r, self.sg, zee),
        }
    }

    /// p/Z(p) at a given pressure.
    fn p_over_z(&self, p: f64) -> f64 {
        p / self.eval_z(p)
    }

    /// Gauss-Legendre integration of 2p/(mu*Z) over [lo, hi] using n-point rule.
    /// Everything stays in Rust — Z-factor, viscosity, quadrature.
    fn gl_integrate(&self, lo: f64, hi: f64, nodes: &[f64], weights: &[f64]) -> f64 {
        let p_mid = (lo + hi) * 0.5;
        let p_half = (hi - lo) * 0.5;
        let mut result = 0.0;

        for (i, &node) in nodes.iter().enumerate() {
            let p_eval = p_mid + p_half * node;
            let zee = self.eval_z(p_eval);
            let ug = self.eval_ug(p_eval, zee);
            let mugz = ug * zee;
            result += weights[i] * 2.0 * p_eval / mugz;
        }

        p_half * result
    }

    /// Solve P/Z -> P for a single target using Newton with bisection fallback.
    fn solve_ponz2p_single(&self, target: f64, rtol: f64) -> Result<f64, String> {
        // Initial guess: p = target (Z ~ 1 at low pressure)
        let mut p = target;
        let p_min = target * 0.1;
        let p_max = target * 5.0;

        // Newton iterations with finite-difference derivative
        let max_newton = 30;
        for _ in 0..max_newton {
            if p < p_min || p > p_max {
                break; // Out of bounds, fall back to bisection
            }
            let f_val = self.p_over_z(p) - target;
            let rel_err = f_val.abs() / target;
            if rel_err < rtol {
                return Ok(p);
            }
            // Central difference for derivative
            let dp = (p * 1e-6).max(0.01);
            let f_plus = self.p_over_z(p + dp) - target;
            let f_minus = self.p_over_z(p - dp) - target;
            let deriv = (f_plus - f_minus) / (2.0 * dp);
            if deriv.abs() < 1e-30 {
                break; // Zero derivative, fall back to bisection
            }
            let step = f_val / deriv;
            p -= step;
        }

        // Bisection fallback
        let mut lo = p_min;
        let mut hi = p_max;
        let f_lo = self.p_over_z(lo) - target;
        let f_hi = self.p_over_z(hi) - target;
        if f_lo * f_hi > 0.0 {
            return Err(two_phase_err(target));
        }

        for _ in 0..100 {
            let mid = (lo + hi) * 0.5;
            let f_mid = self.p_over_z(mid) - target;
            let rel_err = f_mid.abs() / target;
            if rel_err < rtol {
                return Ok(mid);
            }
            if f_hi * f_mid < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Err(two_phase_err(target))
    }
}

fn two_phase_err(target: f64) -> String {
    format!(
        "gas_ponz2p: no single-phase solution exists for P/Z={:.4}. \
         Target may fall in the two-phase region where P/Z is discontinuous.",
        target
    )
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

    // Precompute critical properties — user tc/pc override is method-dependent:
    // DAK/HY+SUT: (tc, pc) replace the mixture Tc/Pc.
    // BNS: (tc, pc) replace only the HC pseudo-component Tc/Pc; inert Tc/Pc
    //      remain at BNS internal constants. LBC params also honor the HC override.
    let user_tc_pc = tc > 0.0 && pc > 0.0;
    let (tpc_sut, ppc_sut) = match method {
        ZMethod::DakSut | ZMethod::HySut => {
            if user_tc_pc {
                (tc, pc)
            } else {
                critical_properties::sutton_wa_internal(sg, co2, h2s, n2)
                    .map_err(pyo3::exceptions::PyValueError::new_err)?
            }
        }
        ZMethod::Bns => (0.0, 0.0), // Not used
    };

    let (tpc_bns, ppc_bns) = match method {
        ZMethod::Bns => {
            if user_tc_pc {
                (tc, pc)
            } else {
                let (t, p, _) = critical_properties::bns_pseudocritical_internal(sg, co2, h2s, n2, h2);
                (t, p)
            }
        }
        _ => (0.0, 0.0),
    };

    // Precompute LBC params for BNS method (propagates HC Tc/Pc override)
    let lbc_p = match method {
        ZMethod::Bns => {
            let (tc_lbc, pc_lbc) = if user_tc_pc { (tc, pc) } else { (0.0, 0.0) };
            Some(gas_viscosity::lbc_params(degf, sg, co2, h2s, n2, h2, tc_lbc, pc_lbc))
        }
        _ => None,
    };

    let state = GasState {
        method, deg_r, sg, tpc_sut, ppc_sut, tpc_bns, ppc_bns,
        co2, h2s, n2, h2, lbc_params: &lbc_p,
    };

    // Two-tier Gauss-Legendre integration
    let result_7 = state.gl_integrate(p1, p2, &GL7_NODES, &GL7_WEIGHTS);
    let result_10 = state.gl_integrate(p1, p2, &GL10_NODES, &GL10_WEIGHTS);

    if result_10.abs() < 1e-30 || (result_10 - result_7).abs() / result_10.abs() < 1e-5 {
        Ok(result_10)
    } else {
        // Split into two subintervals, integrate each with 10-point
        let p_mid = (p1 + p2) * 0.5;
        let r_lo = state.gl_integrate(p1, p_mid, &GL10_NODES, &GL10_WEIGHTS);
        let r_hi = state.gl_integrate(p_mid, p2, &GL10_NODES, &GL10_WEIGHTS);
        Ok(r_lo + r_hi)
    }
}

/// Batch P/Z -> P solver entirely in Rust.
/// Returns a Vec of pressures corresponding to the input P/Z values.
#[pyfunction]
pub fn gas_ponz2p_rust(
    poverz: Vec<f64>,
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
    rtol: f64,
) -> PyResult<Vec<f64>> {
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

    // Precompute critical properties — see gas_dmp_rust for override semantics.
    let user_tc_pc = tc > 0.0 && pc > 0.0;
    let (tpc_sut, ppc_sut) = match method {
        ZMethod::DakSut | ZMethod::HySut => {
            if user_tc_pc {
                (tc, pc)
            } else {
                critical_properties::sutton_wa_internal(sg, co2, h2s, n2)
                    .map_err(pyo3::exceptions::PyValueError::new_err)?
            }
        }
        ZMethod::Bns => (0.0, 0.0),
    };

    let (tpc_bns, ppc_bns) = match method {
        ZMethod::Bns => {
            if user_tc_pc {
                (tc, pc)
            } else {
                let (t, p, _) = critical_properties::bns_pseudocritical_internal(sg, co2, h2s, n2, h2);
                (t, p)
            }
        }
        _ => (0.0, 0.0),
    };

    // P/Z needs no viscosity, so the LBC parameters stay unset here.
    let no_lbc = None;
    let state = GasState {
        method, deg_r, sg, tpc_sut, ppc_sut, tpc_bns, ppc_bns,
        co2, h2s, n2, h2, lbc_params: &no_lbc,
    };

    let mut results = Vec::with_capacity(poverz.len());
    for &target in &poverz {
        let p = state
            .solve_ponz2p_single(target, rtol)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        results.push(p);
    }

    Ok(results)
}
