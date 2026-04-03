/// Oil density and formation volume factor via McCain method.
///
/// Port of oil.py Deno_standing_white_mccainhill (below Pb) and Deno_p_gt_pb
/// (above Pb) plus Bo_mccain.

const TSC: f64 = 60.0; // Standard temperature (deg F)

// ─── Below-Pb density: Standing-Witte-McCain-Hill (1995) ────────────

/// Iterative pseudo-liquid density calculation.
///
/// When `sg_sp > 0`, uses the iterative Eq 3.18b/3.18c approach.
/// When `sg_sp <= 0`, falls back to the apparent-density approach using `sg_g` and `api`.
///
/// Returns oil density in lb/ft3 at pressures at or below bubble point.
pub fn deno_below_pb(
    p: f64,
    degf: f64,
    rs: f64,
    sg_g: f64,
    sg_sp: f64,
    sg_o: f64,
    api: f64,
) -> f64 {
    let rho_po: f64;

    if sg_sp > 0.0 {
        // Iterative path (Eq 3.18b-3.18c)
        let a: [f64; 6] = [-49.8930, 85.0149, -3.70373, 0.0479818, 2.98914, -0.0356888];
        let mut rho_est = (52.8 - 0.01 * rs).max(20.0);
        let mut new_rho = rho_est;
        for _ in 0..100 {
            let rhoa = a[0]
                + a[1] * sg_sp
                + a[2] * sg_sp * rho_est
                + a[3] * sg_sp * rho_est * rho_est
                + a[4] * rho_est
                + a[5] * rho_est * rho_est; // Eq 3.18c
            let denom = 73.71 + rs * sg_sp / rhoa; // Eq 3.18b denominator
            new_rho = (rs * sg_sp + 4600.0 * sg_o) / denom;
            if (rho_est - new_rho).abs() < 1e-8 {
                break;
            }
            rho_est = new_rho;
        }
        rho_po = new_rho;
    } else {
        // Non-iterative fallback using sg_g (Eq 3.17e)
        let rhoa = 38.52 * 10.0_f64.powf(-0.00326 * api)
            + (94.75 - 33.93 * api.log10()) * sg_g.log10();
        rho_po = (rs * sg_g + 4600.0 * sg_o) / (73.71 + rs * sg_g / rhoa);
    }

    // Pressure correction (Eq 3.19d)
    let p_k = p / 1000.0;
    let drho_p = (0.167 + 16.181 * 10.0_f64.powf(-0.0425 * rho_po)) * p_k
        - 0.01 * (0.299 + 263.0 * 10.0_f64.powf(-0.0603 * rho_po)) * p_k * p_k;
    let rho_bs = rho_po + drho_p; // Eq 3.19e

    // Temperature correction (Eq 3.19f)
    let dt = (degf - TSC).max(0.001);
    let drho_t = (0.00302 + 1.505 * rho_bs.powf(-0.951)) * dt.powf(0.938)
        - (0.0216 - 0.0233 * 10.0_f64.powf(-0.0161 * rho_bs)) * dt.powf(0.475);

    rho_bs - drho_t // Eq 3.19g
}

// ─── Above-Pb density: compressibility correction (Eq 3.20) ────────

/// Oil density above bubble point via exponential compressibility correction.
///
/// Uses density at Pb from `deno_below_pb` and the cofb polynomial from
/// McCain Eq 3.13 (compressibility at current pressure).
pub fn deno_above_pb(
    p: f64,
    degf: f64,
    rsb: f64,
    sg_g: f64,
    sg_sp: f64,
    pb: f64,
    sg_o: f64,
    api: f64,
) -> f64 {
    let rho_rb = deno_below_pb(pb, degf, rsb, sg_g, sg_sp, sg_o, api);

    // cofb from McCain Eq 3.13 — polynomial in transformed variables
    let c: [[f64; 6]; 3] = [
        [3.011, -0.0835, 3.51, 0.327, -1.918, 2.52],
        [-2.6254, -0.259, -0.0289, -0.608, -0.642, -2.73],
        [0.497, 0.382, -0.0584, 0.0911, 0.154, 0.429],
    ];
    // Guard against log of non-positive values
    let ln_api = if api > 0.0 { api.ln() } else { 0.0 };
    let ln_sgsp = if sg_sp > 0.0 {
        sg_sp.ln()
    } else if sg_g > 0.0 {
        sg_g.ln()
    } else {
        0.0
    };
    let ln_pb = if pb > 0.0 { pb.ln() } else { 0.0 };
    let p_over_pb = if pb > 0.0 { p / pb } else { 1.0 };
    let ln_p_ratio = if p_over_pb > 0.0 { p_over_pb.ln() } else { 0.0 };
    let ln_rsb = if rsb > 0.0 { rsb.ln() } else { 0.0 };
    let ln_degf = if degf > 0.0 { degf.ln() } else { 0.0 };

    let var: [f64; 6] = [ln_api, ln_sgsp, ln_pb, ln_p_ratio, ln_rsb, ln_degf];

    // Zn[n] = sum over i=0..2 of C[i][n] * var[n]^i
    let mut zp: f64 = 0.0;
    for n in 0..6 {
        let mut zn = 0.0;
        let mut var_pow = 1.0; // var[n]^0 = 1
        for i in 0..3 {
            zn += c[i][n] * var_pow;
            var_pow *= var[n];
        }
        zp += zn;
    }

    let ln_cofb_p = 2.434 + 0.475 * zp + 0.048 * zp * zp - (1.0e6_f64).ln();
    let cofb_p = ln_cofb_p.exp();

    rho_rb * (cofb_p * (p - pb)).exp() // Eq 3.20
}

// ─── Public dispatchers ─────────────────────────────────────────────

/// Oil specific gravity from API.
#[inline]
fn oil_sg(api: f64) -> f64 {
    141.5 / (api + 131.5)
}

/// Oil density via McCain (SWMH) method in lb/ft3.
///
/// Dispatches to `deno_below_pb` when `p <= pb` and `deno_above_pb` otherwise.
/// Either `sg_o` or `api` (or both) must be positive; api takes precedence.
pub fn oil_deno_mccain(
    p: f64,
    degf: f64,
    rs: f64,
    rsb: f64,
    sg_g: f64,
    sg_sp: f64,
    pb: f64,
    sg_o: f64,
    api: f64,
) -> f64 {
    // Resolve sg_o and api consistently (api wins if both supplied)
    let (sg_o_r, api_r) = if api > 0.0 {
        (oil_sg(api), api)
    } else if sg_o > 0.0 {
        (sg_o, 141.5 / sg_o - 131.5)
    } else {
        // Fallback — shouldn't happen with valid inputs
        (0.85, 35.0)
    };

    if p > pb {
        deno_above_pb(p, degf, rsb, sg_g, sg_sp, pb, sg_o_r, api_r)
    } else {
        deno_below_pb(p, degf, rs, sg_g, sg_sp, sg_o_r, api_r)
    }
}

/// Oil formation volume factor via McCain method (rb/STB).
///
/// Bo = (sg_o * 62.372 + 0.01357 * rs * sg_g) / rho_r   (Eq 3.21)
pub fn oil_bo_mccain(
    p: f64,
    degf: f64,
    rs: f64,
    rsb: f64,
    sg_g: f64,
    sg_sp: f64,
    pb: f64,
    sg_o: f64,
    api: f64,
) -> f64 {
    // Resolve sg_o consistently
    let sg_o_r = if api > 0.0 { oil_sg(api) } else { sg_o };
    let sg_g_r = if sg_g > 0.0 { sg_g } else { sg_sp };

    let rho_r = oil_deno_mccain(p, degf, rs, rsb, sg_g, sg_sp, pb, sg_o, api);
    if rho_r <= 0.0 {
        return 1.0; // guard against division by zero
    }
    (sg_o_r * 62.372 + 0.01357 * rs * sg_g_r) / rho_r // Eq 3.21
}

// ─── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: check relative error < tol.
    fn assert_rel(actual: f64, expected: f64, tol: f64, msg: &str) {
        let rel = if expected.abs() > 1e-30 {
            ((actual - expected) / expected).abs()
        } else {
            (actual - expected).abs()
        };
        assert!(
            rel < tol,
            "{}: actual={}, expected={}, rel_err={}",
            msg, actual, expected, rel
        );
    }

    #[test]
    fn test_deno_below_pb_with_sgsp() {
        // Typical black oil: 35 API, sg_sp=0.8, 200 degF, 500 scf/stb, 2000 psia
        let rho = deno_below_pb(2000.0, 200.0, 500.0, 0.8, 0.8, oil_sg(35.0), 35.0);
        // Density should be in reasonable range (30-55 lb/ft3 for live oil)
        assert!(rho > 30.0 && rho < 55.0, "rho={} out of range", rho);
    }

    #[test]
    fn test_deno_below_pb_no_sgsp() {
        // Fallback path with sg_g only
        let rho = deno_below_pb(2000.0, 200.0, 500.0, 0.8, 0.0, oil_sg(35.0), 35.0);
        assert!(rho > 30.0 && rho < 55.0, "rho={} out of range", rho);
    }

    #[test]
    fn test_deno_above_pb() {
        // Above Pb, density should increase with pressure
        let pb = 3000.0;
        let rho_pb = deno_below_pb(pb, 200.0, 600.0, 0.8, 0.8, oil_sg(35.0), 35.0);
        let rho_above = deno_above_pb(5000.0, 200.0, 600.0, 0.8, 0.8, pb, oil_sg(35.0), 35.0);
        assert!(
            rho_above > rho_pb,
            "Above-Pb density ({}) should exceed Pb density ({})",
            rho_above, rho_pb
        );
    }

    #[test]
    fn test_oil_deno_mccain_dispatch() {
        let pb = 3000.0;
        // Below Pb — range check
        let rho_below = oil_deno_mccain(2000.0, 200.0, 400.0, 600.0, 0.8, 0.8, pb, 0.0, 35.0);
        assert!(rho_below > 30.0 && rho_below < 55.0, "rho_below={}", rho_below);

        // At Pb
        let rho_at_pb = oil_deno_mccain(pb, 200.0, 600.0, 600.0, 0.8, 0.8, pb, 0.0, 35.0);
        assert!(rho_at_pb > 30.0 && rho_at_pb < 55.0, "rho_at_pb={}", rho_at_pb);

        // Above Pb — density should exceed density at Pb (compressed liquid)
        let rho_above = oil_deno_mccain(5000.0, 200.0, 600.0, 600.0, 0.8, 0.8, pb, 0.0, 35.0);
        assert!(rho_above > 30.0 && rho_above < 60.0, "rho_above={}", rho_above);
        assert!(rho_above > rho_at_pb,
            "above-Pb density ({}) should exceed density at Pb ({})",
            rho_above, rho_at_pb);
    }

    #[test]
    fn test_oil_bo_mccain_range() {
        let pb = 3000.0;
        // Bo at bubble point — typical range 1.1–1.8
        let bo = oil_bo_mccain(pb, 200.0, 600.0, 600.0, 0.8, 0.8, pb, 0.0, 35.0);
        assert!(bo > 1.0 && bo < 2.5, "bo={} out of range", bo);
    }

    #[test]
    fn test_oil_bo_above_pb_decreases() {
        // Bo should decrease above Pb as oil is compressed
        let pb = 3000.0;
        let bo_pb = oil_bo_mccain(pb, 200.0, 600.0, 600.0, 0.8, 0.8, pb, 0.0, 35.0);
        let bo_hi = oil_bo_mccain(5000.0, 200.0, 600.0, 600.0, 0.8, 0.8, pb, 0.0, 35.0);
        assert!(
            bo_hi < bo_pb,
            "Bo above Pb ({}) should be less than Bo at Pb ({})",
            bo_hi, bo_pb
        );
    }

    #[test]
    fn test_oil_sg_roundtrip() {
        let api = 35.0;
        let sg = oil_sg(api);
        let api_back = 141.5 / sg - 131.5;
        assert_rel(api_back, api, 1e-10, "API roundtrip");
    }
}
