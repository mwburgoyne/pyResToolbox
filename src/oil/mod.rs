use pyo3::prelude::*;

mod density;

/// Oil density via Standing-Witte-McCain-Hill method (lb/ft3).
///
/// Dispatches to below-Pb iterative pseudo-liquid density calculation or
/// above-Pb exponential compressibility correction depending on p vs pb.
///
/// # Arguments
/// * `p`     - Pressure (psia)
/// * `degf`  - Temperature (deg F)
/// * `rs`    - Solution gas ratio at current pressure (scf/stb)
/// * `rsb`   - Solution gas ratio at bubble point (scf/stb)
/// * `sg_g`  - Weighted average surface gas SG (rel. to air)
/// * `sg_sp` - Separator gas SG (rel. to air)
/// * `pb`    - Bubble point pressure (psia)
/// * `sg_o`  - Stock tank oil SG (rel. to water). Ignored when api > 0.
/// * `api`   - Stock tank oil API gravity. Takes precedence over sg_o.
#[pyfunction]
pub fn oil_deno_mccain_rust(
    p: f64,
    degf: f64,
    rs: f64,
    rsb: f64,
    sg_g: f64,
    sg_sp: f64,
    pb: f64,
    sg_o: f64,
    api: f64,
) -> PyResult<f64> {
    Ok(density::oil_deno_mccain(p, degf, rs, rsb, sg_g, sg_sp, pb, sg_o, api))
}

/// Oil formation volume factor via McCain method (rb/STB).
///
/// Bo = (sg_o * 62.372 + 0.01357 * rs * sg_g) / rho_r   (Eq 3.21)
///
/// Uses `oil_deno_mccain` internally for density calculation.
///
/// # Arguments
/// * `p`     - Pressure (psia)
/// * `degf`  - Temperature (deg F)
/// * `rs`    - Solution gas ratio at current pressure (scf/stb)
/// * `rsb`   - Solution gas ratio at bubble point (scf/stb)
/// * `sg_g`  - Weighted average surface gas SG (rel. to air)
/// * `sg_sp` - Separator gas SG (rel. to air)
/// * `pb`    - Bubble point pressure (psia)
/// * `sg_o`  - Stock tank oil SG (rel. to water). Ignored when api > 0.
/// * `api`   - Stock tank oil API gravity. Takes precedence over sg_o.
#[pyfunction]
pub fn oil_bo_mccain_rust(
    p: f64,
    degf: f64,
    rs: f64,
    rsb: f64,
    sg_g: f64,
    sg_sp: f64,
    pb: f64,
    sg_o: f64,
    api: f64,
) -> PyResult<f64> {
    Ok(density::oil_bo_mccain(p, degf, rs, rsb, sg_g, sg_sp, pb, sg_o, api))
}
