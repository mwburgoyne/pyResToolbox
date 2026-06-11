/// Shared VLP segment march for all 4 methods (HB, WG, GRAY, BB) and both
/// fluid types. Direct port of nodal.py _segment_march_gas and
/// _segment_march_oil with per-method gradient callbacks, mirroring the
/// Python scaffold design. nodal.py is the authoritative reference.

use super::constants::*;
use super::friction::serghides_fanning;
use super::holdup_bb::*;
use super::holdup_gray::{gray_effective_roughness, gray_liquid_holdup};
use super::holdup_wg::{wg_friction_gradient_lm, wg_void_fraction};
use super::ift::interfacial_tension;
use super::pvt_helpers::*;
use super::static_column::{static_gas_column_pressure, static_oil_column_pressure};

/// Divergence guard message - must match the Python RuntimeError text in
/// nodal.py exactly.
pub const DIVERGED_MSG: &str = "VLP march diverged: pressure fell below atmospheric - the \
specified rates are not physically sustainable for this geometry";

/// Per-step PVT and flow state passed to the gradient callbacks.
/// Mirrors the state dict built by the Python segment marches; only the
/// fields consumed by the four gradient functions are carried.
pub struct SegmentState {
    pub p_avg: f64,    // Average segment pressure (psia)
    pub mu_g: f64,     // Gas viscosity (cP)
    pub rho_g: f64,    // Gas density (lb/ft^3)
    pub rho_l: f64,    // Liquid density (lb/ft^3)
    pub rho_ns: f64,   // No-slip mixture density (lb/ft^3)
    pub v_sg: f64,     // Superficial gas velocity (ft/s)
    pub v_sl: f64,     // Superficial liquid velocity (ft/s)
    pub v_m: f64,      // Mixture velocity (ft/s)
    pub lambda_l: f64, // No-slip liquid holdup
    pub sigma: f64,    // Interfacial tension (dyne/cm)
    pub mu_l: f64,     // Liquid viscosity (cP)
    pub diam_ft: f64,  // Pipe internal diameter (ft)
    pub rough_ft: f64, // Pipe roughness (ft)
    pub area: f64,     // Pipe cross-sectional area (ft^2)
    pub injection: bool,
    pub theta: f64, // Angle from horizontal (radians)
    pub mflow_g: f64, // Gas mass flow (lbm/s)
    pub mflow_l: f64, // Liquid mass flow (lbm/s)
    pub mflow_total: f64, // Total mass flow (lbm/s)
    pub ql_loc: f64,  // Local liquid rate (STB/d)
    pub lsg_loc: f64, // Local liquid specific gravity
    pub tid: f64,     // Pipe internal diameter (inches)
    pub rough: f64,   // Pipe roughness (inches)
}

pub type GradientFn = fn(&SegmentState) -> f64;

/// Number of march segments for a given length (Python _calc_segments,
/// min_seg_ft = 100).
fn calc_segments(length: f64) -> usize {
    let mut ndiv = 100usize;
    let seg_len = length / ndiv as f64;
    if seg_len < 100.0 {
        ndiv = (length / 100.0) as usize + 1;
        if ndiv < 1 {
            ndiv = 1;
        }
    }
    ndiv
}

/// Linear CGR interpolation for condensate dropout (Python _condensate_dropout).
fn condensate_dropout(
    cgr: f64, qg_mmscfd: f64, p_avg: f64, pr: f64,
    osg: f64, qw_bwpd: f64, wsg: f64,
) -> (f64, f64, f64, f64) {
    let cgr_local = if pr > 14.7 {
        cgr * (pr - p_avg).max(0.0) / (pr - 14.7)
    } else {
        cgr
    };
    let qo_local = cgr_local * qg_mmscfd;
    let ql_local = qo_local + qw_bwpd;
    let lsg_local = if ql_local > 0.0 {
        (qo_local * osg + qw_bwpd * wsg) / ql_local
    } else {
        osg // No liquid; value won't be used meaningfully
    };
    (cgr_local, qo_local, ql_local, lsg_local)
}

/// Condensate viscosity (Python _condensate_vis).
fn condensate_vis(
    pr: f64, cgr_local: f64, gsg: f64, api: f64, temp_f: f64,
    p_avg: f64, oil_vis: f64,
) -> f64 {
    if pr > 14.7 && cgr_local > 0.01 {
        oil_viscosity_full(gsg, api, temp_f, 0.0, 14.7, p_avg, 1.0, 1.0)
    } else {
        oil_vis
    }
}

// ============================================================================
//  Gradient callbacks (psi/ft), one per VLP method
// ============================================================================

/// Hagedorn-Brown gradient: polynomial holdup, Orkiszewski bubble flow,
/// Serghides friction. Port of nodal.py _hb_gradient_gas.
fn hb_gradient(s: &SegmentState) -> f64 {
    let ul = s.v_sl;
    let ugas = s.v_sg;
    let ift = s.sigma;
    let rho_l = s.rho_l;

    let nvl = DG_VEL * ul * (rho_l / ift).powf(0.25);
    let nvg = DG_VEL * ugas * (rho_l / ift).powf(0.25);
    let nd = DG_DIAM * s.diam_ft * (rho_l / ift).powf(0.5);

    let mul = if s.ql_loc > 0.0 { s.mu_l } else { 1.0 };

    let mut nl = DG_VISC * mul * (1.0 / (rho_l * ift.powi(3))).powf(0.25);
    if nl <= 0.0 {
        nl = 1e-8;
    }

    let x1 = log10_safe(nl) + 3.0;
    let cnl = 10.0_f64.powf(
        CNL_C0 + CNL_C1 * x1 + CNL_C2 * x1 * x1 + CNL_C3 * x1.powi(3) + CNL_C4 * x1.powi(4),
    );

    let nvg_safe = nvg.max(1e-10);
    let f1 = nvl * s.p_avg.powf(0.1) * cnl / (nvg_safe.powf(F1_VG_EXP) * F1_P_ATM.powf(0.1) * nd);

    let lf1 = log10_safe(f1) + 6.0;
    let mut ylonsi =
        YL_C0 + YL_C1 * lf1 + YL_C2 * lf1.powi(2) + YL_C3 * lf1.powi(3) + YL_C4 * lf1.powi(4);
    ylonsi = ylonsi.max(0.0);

    let mut f2 = nvg * nl.powf(F2_VISC_EXP) / nd.powf(F2_DIAM_EXP);
    f2 = f2.max(SC_F2_CLAMP);

    let mut si = SC_C0 + SC_C1 * f2 + SC_C2 * f2 * f2 + SC_C3 * f2.powi(3) + SC_C4 * f2.powi(4);
    if f2 <= SC_F2_FLOOR {
        si = 1.0;
    }

    let mut yl = clamp(si * ylonsi, 0.0, 1.0);

    // Minimum holdup from mass fraction
    let rho_g_sg = s.rho_g / 62.37;
    let mflow_lpd = s.mflow_l * SEC_PER_DAY;
    let mflow_gpd = s.mflow_g * SEC_PER_DAY;
    let mflow_total_pd = mflow_lpd + mflow_gpd;
    let mass_frac_liq = if mflow_total_pd > 0.0 {
        mflow_lpd / mflow_total_pd
    } else {
        0.0
    };
    let min_l = if (rho_g_sg + s.lsg_loc) > 0.0 {
        mass_frac_liq * rho_g_sg / (rho_g_sg + s.lsg_loc)
    } else {
        0.0
    };
    if yl < min_l {
        yl = min_l;
    }

    // Orkiszewski bubble flow correction
    let vm = ugas + ul;
    let vs = ORK_VS;
    let mut lb = ORK_LB_A + ORK_LB_B * vm * vm / s.diam_ft;
    if lb < ORK_LB_MIN {
        lb = ORK_LB_MIN;
    }
    if vm > 0.0 {
        let b_ratio = ugas / vm;
        if b_ratio < lb {
            let disc = (1.0 + vm / vs).powi(2) - 4.0 * ugas / vs;
            if disc >= 0.0 {
                yl = 1.0 - 0.5 * (1.0 + vm / vs - disc.sqrt());
                yl = clamp(yl, 0.0, 1.0);
            }
        }
    }

    // Reynolds and friction (mass-flow basis)
    let mflow_pd = mflow_total_pd;
    let mut nre =
        HB_RE_K * (mflow_pd / SEC_PER_DAY) / (s.diam_ft * mul.powf(yl) * s.mu_g.powf(1.0 - yl));
    if nre < 2100.0 {
        nre = 2100.0;
    }

    let eond = s.rough / s.tid;
    let f = serghides_fanning(nre, eond);

    let rho_avg = yl * rho_l + (1.0 - yl) * s.rho_g;

    let fric_term = f * mflow_pd.powi(2) / HB_FRIC_K / s.diam_ft.powi(5) / rho_avg;
    let sign = if s.injection { -1.0 } else { 1.0 };
    (rho_avg * s.theta.sin() + sign * fric_term) / IN2_PER_FT2
}

/// Woldesemayat-Ghajar gradient: drift-flux void fraction, Chisholm friction
/// (SI). Port of nodal.py _wg_gradient_gas.
fn wg_gradient(s: &SegmentState) -> f64 {
    // Convert oilfield state to SI
    let diam_m = s.tid * IN_TO_M;
    let rough_m = s.rough * IN_TO_M;
    let area_m2 = std::f64::consts::PI * diam_m.powi(2) / 4.0;
    let rho_g_kgm3 = s.rho_g * LBFT3_TO_KGM3;
    let rho_l_kgm3 = s.rho_l * LBFT3_TO_KGM3;
    let mflow_g_kgs = s.mflow_g * LB_TO_KG;
    let mflow_l_kgs = s.mflow_l * LB_TO_KG;
    let u_sg = mflow_g_kgs / (rho_g_kgm3 * area_m2);
    let u_sl = mflow_l_kgs / (rho_l_kgm3 * area_m2);
    let sigma_nm = s.sigma * DYNECM_TO_NM;
    let p_sys_pa = s.p_avg * PSI_TO_PA;

    let alpha_g = wg_void_fraction(
        u_sg, u_sl, rho_g_kgm3, rho_l_kgm3, sigma_nm, diam_m, s.theta, p_sys_pa,
    );
    let rho_mix = alpha_g * rho_g_kgm3 + (1.0 - alpha_g) * rho_l_kgm3;
    let dpdz_hydro = rho_mix * G_SI * s.theta.sin();

    let mu_g_pas = s.mu_g * CP_TO_PAS;
    let mu_l_pas = s.mu_l * CP_TO_PAS;
    let dpdz_fric = wg_friction_gradient_lm(
        mflow_g_kgs, mflow_l_kgs, rho_g_kgm3, rho_l_kgm3, mu_g_pas, mu_l_pas, diam_m, rough_m,
    );

    let sign = if s.injection { -1.0 } else { 1.0 };
    let dpdz_pam = dpdz_hydro + sign * dpdz_fric;
    dpdz_pam / (PSI_TO_PA / FT_TO_M)
}

/// Gray method gradient: holdup via effective roughness, acceleration term.
/// Port of nodal.py _gray_gradient_gas.
fn gray_gradient(s: &SegmentState) -> f64 {
    let hl = gray_liquid_holdup(s.v_sl, s.v_sg, s.rho_l, s.rho_g, s.sigma, s.diam_ft, s.lambda_l);
    let rho_s = s.rho_l * hl + s.rho_g * (1.0 - hl);

    let eps_eff = gray_effective_roughness(s.rough_ft, s.sigma, s.rho_ns, s.v_sl, s.v_sg);
    let eps_d = eps_eff / s.diam_ft;

    let mu_ns = s.mu_l * s.lambda_l + s.mu_g * (1.0 - s.lambda_l);
    let mu_ns_lbfts = mu_ns * CP_TO_LBFTS;
    let n_re = if mu_ns_lbfts > 0.0 {
        s.rho_ns * s.v_m * s.diam_ft / mu_ns_lbfts
    } else {
        0.0
    };

    let f_moody = 4.0 * serghides_fanning(n_re, eps_d);

    let gm = s.mflow_total / s.area;
    let dpdz_hydro = rho_s * s.theta.sin() / IN2_PER_FT2;
    let dpdz_fric = if s.rho_ns > 0.0 {
        f_moody * gm.powi(2) / (2.0 * GC * s.diam_ft * s.rho_ns * IN2_PER_FT2)
    } else {
        0.0
    };
    let ek = gm * s.v_sg / (GC * s.p_avg * IN2_PER_FT2);
    let denom = (1.0 - ek).max(0.1);

    let sign = if s.injection { -1.0 } else { 1.0 };
    (dpdz_hydro + sign * dpdz_fric) / denom
}

/// Beggs and Brill gradient: flow pattern map, Payne correction, two-phase
/// friction. Port of nodal.py _bb_gradient_gas.
fn bb_gradient(s: &SegmentState) -> f64 {
    let froude = s.v_m.powi(2) / (G_FT * s.diam_ft);
    let (pattern, trans_a) = bb_flow_pattern(froude, s.lambda_l);

    let mut hl0 = if pattern == BB_TRANSITION {
        let hl0_seg = bb_horizontal_holdup(s.lambda_l, froude, BB_SEGREGATED);
        let hl0_int = bb_horizontal_holdup(s.lambda_l, froude, BB_INTERMITTENT);
        trans_a * hl0_seg + (1.0 - trans_a) * hl0_int
    } else {
        bb_horizontal_holdup(s.lambda_l, froude, pattern)
    };

    // Payne et al. (1979), JPT 31(9): uphill liquid holdup correction factor
    // 0.924, applied to all flow patterns (holdup floored at no-slip lambda_l)
    hl0 *= BB_PAYNE;
    hl0 = hl0.max(s.lambda_l);

    // Liquid velocity number NLV = 1.938 * vsl * (rho_l/sigma)^0.25 with sigma
    // in dyne/cm (Beggs and Brill 1973; same form and units as Hagedorn-Brown)
    let n_lv = if s.sigma > 0.0 {
        DG_VEL * s.v_sl * (s.rho_l / s.sigma).powf(0.25)
    } else {
        0.0
    };

    let hl_theta = if pattern == BB_TRANSITION {
        let hl_seg = bb_inclination_correction(
            bb_horizontal_holdup(s.lambda_l, froude, BB_SEGREGATED) * BB_PAYNE,
            s.lambda_l, n_lv, froude, BB_SEGREGATED, s.theta,
        );
        let hl_int = bb_inclination_correction(
            bb_horizontal_holdup(s.lambda_l, froude, BB_INTERMITTENT) * BB_PAYNE,
            s.lambda_l, n_lv, froude, BB_INTERMITTENT, s.theta,
        );
        trans_a * hl_seg + (1.0 - trans_a) * hl_int
    } else {
        bb_inclination_correction(hl0, s.lambda_l, n_lv, froude, pattern, s.theta)
    };

    let hl_theta = clamp(hl_theta, s.lambda_l, 1.0);
    let rho_s = s.rho_l * hl_theta + s.rho_g * (1.0 - hl_theta);

    let mu_ns = s.mu_l * s.lambda_l + s.mu_g * (1.0 - s.lambda_l);
    let mu_ns_lbfts = mu_ns * CP_TO_LBFTS;
    let n_re = if mu_ns_lbfts > 0.0 {
        s.rho_ns * s.v_m * s.diam_ft / mu_ns_lbfts
    } else {
        0.0
    };

    let eps_d = s.rough_ft / s.diam_ft;
    let f_ns = serghides_fanning(n_re, eps_d);
    let f_tp = bb_two_phase_friction(f_ns, s.lambda_l, hl_theta);

    let dpdz_hydro = rho_s * s.theta.sin() / IN2_PER_FT2;
    let dpdz_fric = if s.rho_ns > 0.0 {
        4.0 * f_tp * s.rho_ns * s.v_m.powi(2) / (2.0 * GC * s.diam_ft * IN2_PER_FT2)
    } else {
        0.0
    };
    let ek = rho_s * s.v_m * s.v_sg / (GC * s.p_avg * IN2_PER_FT2);
    let denom = (1.0 - ek).max(0.1);

    let sign = if s.injection { -1.0 } else { 1.0 };
    (dpdz_hydro + sign * dpdz_fric) / denom
}

// ============================================================================
//  Shared gas segment march (Python _segment_march_gas)
// ============================================================================

#[allow(clippy::too_many_arguments)]
fn segment_march_gas(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64, gradient: GradientFn,
) -> Result<f64, String> {
    let (tc, pc) = sutton_tc_pc(gsg);
    let osg = 141.5 / (api + 131.5);
    let total_mass = RHO_AIR_STC * gsg * qg_mmscfd * 1e6
        + osg * RHO_FW * cgr * qg_mmscfd * FT3_PER_BBL
        + wsg * RHO_FW * qw_bwpd * FT3_PER_BBL;
    if total_mass < 1e-6 || qg_mmscfd < 0.001 {
        return Ok(static_gas_column_pressure(thp, length, tht, bht, gsg, theta));
    }

    let diam_ft = tid / 12.0;
    let rough_ft = rough / 12.0;
    let area = std::f64::consts::PI * diam_ft.powi(2) / 4.0;

    let ndiv = calc_segments(length);
    let seg_len = length / ndiv as f64;

    let mflow_g = RHO_AIR_STC * gsg * qg_mmscfd * 1e6 / SEC_PER_DAY;
    let mflow_w = wsg * RHO_FW * qw_bwpd * FT3_PER_BBL / SEC_PER_DAY;

    let mut p_psia = thp;

    for i in 1..=ndiv {
        let frac = (i as f64 - 0.5) / ndiv as f64;
        let temp_f = tht + (bht - tht) * frac;
        let temp_r = temp_f + 459.67;
        let mut p_est = p_psia;

        for _ in 0..2 {
            let p_avg = ((p_psia + p_est) / 2.0).max(14.7);

            let (cgr_loc, qo_loc, ql_loc, lsg_loc) =
                condensate_dropout(cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg);

            let mflow_o = osg * RHO_FW * qo_loc * FT3_PER_BBL / SEC_PER_DAY;
            let mflow_l = mflow_o + mflow_w;
            let mflow_total = mflow_g + mflow_l;

            let oil_vis_loc = condensate_vis(pr, cgr_loc, gsg, api, temp_f, p_avg, oil_vis);

            let zee = z_factor(gsg, temp_f, p_avg, tc, pc);
            let mu_g = gas_viscosity(gsg, temp_f, p_avg, zee, tc, pc);

            let rho_g = MW_AIR * gsg * p_avg / (zee * R_GAS * temp_r);
            let rho_l = lsg_loc * RHO_FW;

            let v_sg = (mflow_g / rho_g.max(1e-10)) / area;
            let v_sl = (mflow_l / rho_l.max(1e-10)) / area;
            let v_m = v_sg + v_sl;
            let lambda_l = if v_m > 1e-10 { v_sl / v_m } else { 0.0 };
            let rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l);

            let water_visc = water_viscosity(p_avg, temp_f, 0.0);
            let mu_l = if ql_loc > 0.0 {
                (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc
            } else {
                oil_vis_loc
            };

            let rs_est = standing_rs(gsg, p_avg, temp_f, api);
            let sigma = interfacial_tension(
                p_avg, temp_f, api, rs_est,
                mflow_o * LB_TO_KG, mflow_w * LB_TO_KG, mflow_g * LB_TO_KG,
            );

            let s = SegmentState {
                p_avg, mu_g, rho_g, rho_l, rho_ns, v_sg, v_sl, v_m, lambda_l,
                sigma, mu_l, diam_ft, rough_ft, area, injection, theta,
                mflow_g, mflow_l, mflow_total, ql_loc, lsg_loc, tid, rough,
            };

            let dpdz = gradient(&s);
            p_est = p_psia + dpdz * seg_len;
        }

        p_psia = p_est;
        if p_psia < 1.0 {
            return Err(DIVERGED_MSG.to_string());
        }
    }

    Ok(p_psia)
}

// ============================================================================
//  Shared oil segment march (Python _segment_march_oil)
// ============================================================================

#[allow(clippy::too_many_arguments)]
fn segment_march_oil(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qt_stbpd: f64, gor: f64, wc: f64, pb: f64, rsb: f64, sgsp: f64,
    rsb_scale: f64, injection: bool, theta: f64, gradient: GradientFn,
    vis_frac: f64, rsb_frac: f64,
) -> Result<f64, String> {
    let (tc, pc) = sutton_tc_pc(gsg);
    let wc_adj = wc.max(1e-9);
    if qt_stbpd < 1e-7 {
        return Ok(static_oil_column_pressure(
            thp, length, tht, bht, wc_adj, wsg, api, sgsp, pb, rsb, rsb_scale, theta,
        ));
    }

    let qo = qt_stbpd * (1.0 - wc);
    let qw = qt_stbpd * wc;
    let osg = 141.5 / (api + 131.5);
    let rsb_for_calc = rsb / rsb_scale;

    let diam_ft = tid / 12.0;
    let rough_ft = rough / 12.0;
    let area = std::f64::consts::PI * diam_ft.powi(2) / 4.0;

    let ndiv = calc_segments(length);
    let seg_len = length / ndiv as f64;

    let mut p_psia = thp;

    for i in 1..=ndiv {
        let frac = (i as f64 - 0.5) / ndiv as f64;
        let temp_f = tht + (bht - tht) * frac;
        let temp_r = temp_f + 459.67;
        let mut p_est = p_psia;

        for _ in 0..2 {
            let p_avg = ((p_psia + p_est) / 2.0).max(14.7);

            let rs_local = velarde_rs(sgsp, api, temp_f, pb, rsb_for_calc, p_avg) * rsb_scale;
            let free_gas = (gor - rs_local).max(0.0);
            let qg_mmscfd = (free_gas * qo / 1e6).max(1e-9);

            let oil_vis_seg =
                oil_viscosity_full(sgsp, api, temp_f, rsb, pb, p_avg, vis_frac, rsb_frac);
            let rho_oil = oil_density_mccain(rs_local, sgsp, osg, p_avg.min(pb), temp_f);

            let zee = z_factor(gsg, temp_f, p_avg, tc, pc);
            let mu_g = gas_viscosity(gsg, temp_f, p_avg, zee, tc, pc);

            let rho_g = MW_AIR * gsg * p_avg / (zee * R_GAS * temp_r);

            let mflow_o = osg * RHO_FW * qo * FT3_PER_BBL / SEC_PER_DAY;
            let mflow_w = wsg * RHO_FW * qw * FT3_PER_BBL / SEC_PER_DAY;
            let mflow_g = RHO_AIR_STC * gsg * qg_mmscfd * 1e6 / SEC_PER_DAY;
            let mflow_l = mflow_o + mflow_w;
            let mflow_total = mflow_g + mflow_l;

            let ql = qo + qw;
            let rho_w = wsg * RHO_FW;
            // Liquid mixture density uses the live-oil (McCain) density for the
            // oil fraction - matches Python _segment_march_oil
            let rho_l = if ql > 0.0 {
                (qo * rho_oil + qw * rho_w) / ql
            } else {
                rho_oil
            };
            let lsg = if ql > 0.0 {
                (qo * osg + qw * wsg) / ql
            } else {
                osg
            };

            let v_sg = (mflow_g / rho_g.max(1e-10)) / area;
            let v_sl = (mflow_l / rho_l.max(1e-10)) / area;
            let v_m = v_sg + v_sl;
            let lambda_l = if v_m > 1e-10 { v_sl / v_m } else { 0.0 };
            let rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l);

            let water_visc = water_viscosity(p_avg, temp_f, 0.0);
            let mu_l = if ql > 0.0 {
                (qo * oil_vis_seg + qw * water_visc) / ql
            } else {
                oil_vis_seg
            };

            let sigma = interfacial_tension(
                p_avg, temp_f, api, rs_local,
                mflow_o * LB_TO_KG, mflow_w * LB_TO_KG, mflow_g * LB_TO_KG,
            );

            let s = SegmentState {
                p_avg, mu_g, rho_g, rho_l, rho_ns, v_sg, v_sl, v_m, lambda_l,
                sigma, mu_l, diam_ft, rough_ft, area, injection, theta,
                mflow_g, mflow_l, mflow_total,
                ql_loc: ql, lsg_loc: lsg, tid, rough,
            };

            let dpdz = gradient(&s);
            p_est = p_psia + dpdz * seg_len;
        }

        p_psia = p_est;
        if p_psia < 1.0 {
            return Err(DIVERGED_MSG.to_string());
        }
    }

    Ok(p_psia)
}

// ============================================================================
//  Public entry points (4 methods x 2 fluid types)
// ============================================================================

macro_rules! gas_entry {
    ($name:ident, $gradient:expr) => {
        #[allow(clippy::too_many_arguments)]
        pub fn $name(
            thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
            length: f64, tht: f64, bht: f64, wsg: f64,
            qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
            injection: bool, pr: f64, theta: f64,
        ) -> Result<f64, String> {
            segment_march_gas(
                thp, api, gsg, tid, rough, length, tht, bht, wsg,
                qg_mmscfd, cgr, qw_bwpd, oil_vis, injection, pr, theta, $gradient,
            )
        }
    };
}

macro_rules! oil_entry {
    ($name:ident, $gradient:expr) => {
        #[allow(clippy::too_many_arguments)]
        pub fn $name(
            thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
            length: f64, tht: f64, bht: f64, wsg: f64,
            qt_stbpd: f64, gor: f64, wc: f64, pb: f64, rsb: f64, sgsp: f64,
            rsb_scale: f64, injection: bool, theta: f64,
            vis_frac: f64, rsb_frac: f64,
        ) -> Result<f64, String> {
            segment_march_oil(
                thp, api, gsg, tid, rough, length, tht, bht, wsg,
                qt_stbpd, gor, wc, pb, rsb, sgsp, rsb_scale, injection, theta,
                $gradient, vis_frac, rsb_frac,
            )
        }
    };
}

gas_entry!(hb_fbhp_gas, hb_gradient);
gas_entry!(wg_fbhp_gas, wg_gradient);
gas_entry!(gray_fbhp_gas, gray_gradient);
gas_entry!(bb_fbhp_gas, bb_gradient);

oil_entry!(hb_fbhp_oil, hb_gradient);
oil_entry!(wg_fbhp_oil, wg_gradient);
oil_entry!(gray_fbhp_oil, gray_gradient);
oil_entry!(bb_fbhp_oil, bb_gradient);
