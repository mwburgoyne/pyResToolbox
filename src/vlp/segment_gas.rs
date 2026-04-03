/// Gas well VLP segment loops for all 4 methods.
/// Direct port of nodal.py _hb_fbhp_gas, _wg_fbhp_gas, _gray_fbhp_gas, _bb_core_gas.

use super::pvt_helpers::*;
use super::ift::*;
use super::friction::serghides_fanning;
use super::holdup_hb::*;
use super::holdup_wg::*;
use super::holdup_gray::*;
use super::holdup_bb::*;
use super::static_column::static_gas_column_pressure;

const MW_AIR: f64 = 28.97;
const G_FT: f64 = 32.174;
const GC: f64 = 32.174;
const G_SI: f64 = 9.80665;
const PSI_TO_PA: f64 = 6894.757;
const FT_TO_M: f64 = 0.3048;
const IN_TO_M: f64 = 0.0254;
const LBFT3_TO_KGM3: f64 = 16.01846;
const CP_TO_PAS: f64 = 0.001;
const DYNECM_TO_NM: f64 = 0.001;

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
        osg
    };
    (cgr_local, qo_local, ql_local, lsg_local)
}

fn condensate_vis(
    pr: f64, cgr_local: f64, gsg: f64, api: f64, temp_f: f64,
    p_avg: f64, default_vis: f64,
) -> f64 {
    if pr > 14.7 && cgr_local > 0.01 {
        oil_viscosity_full(gsg, api, temp_f, 0.0, 14.7, p_avg, 1.0, 1.0)
    } else {
        default_vis
    }
}

// ============================================================================
//  Hagedorn-Brown Gas
// ============================================================================
pub fn hb_fbhp_gas(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> f64 {
    let (tc, pc) = sutton_tc_pc(gsg);
    let osg = 141.5 / (api + 131.5);
    let total_mass = 0.0765 * gsg * qg_mmscfd * 1e6
        + osg * 62.4 * cgr * qg_mmscfd * 5.615
        + wsg * 62.4 * qw_bwpd * 5.615;
    if total_mass < 1e-6 || qg_mmscfd < 0.05 {
        return static_gas_column_pressure(thp, length, tht, bht, gsg, theta);
    }

    let area = std::f64::consts::PI * tid * tid / 4.0 / 144.0;
    let ndiv = calc_segments(length);
    let seg_len = length / ndiv as f64;

    // Build temperature array
    let mut depth = vec![0.0; ndiv + 1];
    let mut dz = vec![0.0; ndiv + 1];
    let mut incr = vec![0.0; ndiv + 1];
    let mut temp = vec![0.0; ndiv + 1];
    let mut p = vec![0.0; ndiv + 1];

    temp[0] = tht;
    let mut midpoint = 0.0;

    for i in 0..=ndiv {
        depth[i] = i as f64 * seg_len;
        if i == 0 {
            incr[i] = seg_len;
        } else {
            incr[i] = depth[i] - depth[i - 1];
            dz[i] = (incr[i - 1] + incr[i]) / 2.0;
        }
        if i == ndiv {
            incr[i] = 0.5 * seg_len;
            dz[i] = incr[i - 1] + incr[i];
        }
        midpoint += dz[i];
        temp[i] = (bht - tht) * midpoint / length + tht + 459.67;
    }

    p[0] = thp;

    for i in 1..=ndiv {
        let temp_f_i = temp[i] - 459.67;
        let mut p_est = p[i - 1];

        for _ in 0..2 {
            let p_avg = ((p[i - 1] + p_est) / 2.0).max(14.7);

            let (cgr_loc, qo_loc, ql_loc, lsg_loc) =
                condensate_dropout(cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg);
            let ul = ql_loc * 5.615 / 86400.0 / area;

            let oil_vis_loc =
                condensate_vis(pr, cgr_loc, gsg, api, temp_f_i, p_avg, oil_vis);

            let zee = z_factor(gsg, temp_f_i, p_avg, tc, pc);

            let mflow_o = osg * 62.4 * qo_loc * 5.615;
            let mflow_w = wsg * 62.4 * qw_bwpd * 5.615;
            let mflow_g = 0.0765 * gsg * qg_mmscfd * 1e6;

            let water_visc = water_viscosity(p_avg, temp_f_i, 0.0);

            let rs_est = standing_rs(gsg, p_avg, temp_f_i, api);
            let ift_val = interfacial_tension(
                p_avg, temp_f_i, api, rs_est, mflow_o, mflow_w, mflow_g,
            );

            let ugas = qg_mmscfd * 1e6
                / (p_avg * 35.3741 / (zee * temp[i]))
                / 86400.0
                / area;

            let mul = if ql_loc > 0.0 {
                (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / (qo_loc + qw_bwpd)
            } else {
                1.0
            };

            let (mut yl, _nvl, _nvg, _nd) =
                hb_holdup(ul, ugas, tid, lsg_loc, ift_val, mul, p_avg);

            // Minimum holdup from mass fraction
            let rho = 29.0 * gsg * p_avg / (temp[i] * zee * 10.732 * 62.37);
            let mflow = lsg_loc * 62.4 * ql_loc * 5.615 + 0.0765 * gsg * qg_mmscfd * 1e6;
            let mass_frac_liq = if mflow > 0.0 {
                (lsg_loc * 62.4 * ql_loc * 5.615) / mflow
            } else {
                0.0
            };
            let min_l = mass_frac_liq * rho * 62.37 / (rho * 62.37 + lsg_loc * 62.37);
            if yl < min_l {
                yl = min_l;
            }

            // Orkiszewski correction
            yl = orkiszewski_correction(yl, ugas, ul, tid);

            let mug = gas_viscosity(gsg, temp_f_i, p_avg, zee, tc, pc);
            let mut nre = 96778.0 * (mflow / 86400.0)
                / ((tid / 12.0) * mul.powf(yl) * mug.powf(1.0 - yl));
            if nre < 2100.0 {
                nre = 2100.0;
            }

            let eond = rough / tid;
            let f = serghides_fanning(nre, eond);

            let rho_g = MW_AIR * gsg * p_avg / (zee * 10.73 * temp[i]);
            let rho_avg = yl * lsg_loc * 62.4 + (1.0 - yl) * rho_g;

            let fric_term =
                f * mflow * mflow / 7.413e10 / (tid / 12.0).powi(5) / rho_avg;
            let dpdz = (rho_avg * theta.sin()
                + if injection { -fric_term } else { fric_term })
                / 144.0;

            p_est = p[i - 1] + dpdz * incr[i];
        }

        p[i] = p_est;
    }

    p[ndiv]
}

// ============================================================================
//  Woldesemayat-Ghajar Gas
// ============================================================================
pub fn wg_fbhp_gas(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> f64 {
    let (tc, pc) = sutton_tc_pc(gsg);
    let osg = 141.5 / (api + 131.5);
    let total_mass = 0.0765 * gsg * qg_mmscfd * 1e6
        + osg * 62.4 * cgr * qg_mmscfd * 5.615
        + wsg * 62.4 * qw_bwpd * 5.615;
    if total_mass < 1e-6 || qg_mmscfd < 0.05 {
        return static_gas_column_pressure(thp, length, tht, bht, gsg, theta);
    }

    let diam_m = tid * IN_TO_M;
    let rough_m = rough * IN_TO_M;
    let area_m2 = std::f64::consts::PI * diam_m * diam_m / 4.0;

    let m_flow_g_kgs = 0.0765 * gsg * qg_mmscfd * 1e6 * 0.453592 / 86400.0;
    let m_flow_w_kgs = wsg * 62.4 * qw_bwpd * 5.615 * 0.453592 / 86400.0;

    let ndiv = calc_segments(length);
    let seg_len_ft = length / ndiv as f64;

    let mut p_psia = thp;

    for i in 1..=ndiv {
        let frac = (i as f64 - 0.5) / ndiv as f64;
        let temp_f_i = tht + (bht - tht) * frac;
        let mut p_est = p_psia;

        for _ in 0..2 {
            let mut p_avg = (p_psia + p_est) / 2.0;
            if p_avg < 14.7 {
                p_avg = 14.7;
            }

            let (cgr_loc, qo_loc, ql_loc, lsg_loc) =
                condensate_dropout(cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg);

            let m_flow_o_kgs = osg * 62.4 * qo_loc * 5.615 * 0.453592 / 86400.0;
            let m_flow_l_kgs = m_flow_o_kgs + m_flow_w_kgs;
            let rho_l_kgm3 = lsg_loc * 62.4 * LBFT3_TO_KGM3;
            let u_sl = m_flow_l_kgs / (rho_l_kgm3 * area_m2);

            let oil_vis_loc =
                condensate_vis(pr, cgr_loc, gsg, api, temp_f_i, p_avg, oil_vis);

            let zee = z_factor(gsg, temp_f_i, p_avg, tc, pc);
            let mu_g_cp = gas_viscosity(gsg, temp_f_i, p_avg, zee, tc, pc);

            let temp_r = temp_f_i + 459.67;
            let rho_g_lbft3 = MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r);
            let rho_g_kgm3 = rho_g_lbft3 * LBFT3_TO_KGM3;

            let qg_actual_m3s = m_flow_g_kgs / rho_g_kgm3;
            let u_sg = qg_actual_m3s / area_m2;

            let water_visc = water_viscosity(p_avg, temp_f_i, 0.0);
            let mu_l = if ql_loc > 0.0 {
                (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc * CP_TO_PAS
            } else {
                oil_vis_loc * CP_TO_PAS
            };

            let rs_est = standing_rs(gsg, p_avg, temp_f_i, api);
            let ift_val = interfacial_tension(
                p_avg, temp_f_i, api, rs_est,
                m_flow_o_kgs, m_flow_w_kgs, m_flow_g_kgs,
            );
            let sigma_nm = ift_val * DYNECM_TO_NM;
            let p_sys_pa = p_avg * PSI_TO_PA;

            let alpha_g = wg_void_fraction(
                u_sg, u_sl, rho_g_kgm3, rho_l_kgm3,
                sigma_nm, diam_m, theta, p_sys_pa,
            );
            let liq_holdup = 1.0 - alpha_g;

            let rho_mix = alpha_g * rho_g_kgm3 + liq_holdup * rho_l_kgm3;
            let dpdz_hydro = rho_mix * G_SI * theta.sin();

            let mu_g_pas = mu_g_cp * CP_TO_PAS;
            let dpdz_fric = wg_friction_gradient_lm(
                m_flow_g_kgs, m_flow_l_kgs, rho_g_kgm3, rho_l_kgm3,
                mu_g_pas, mu_l, diam_m, rough_m,
            );

            let dpdz_total_pam =
                dpdz_hydro + if injection { -dpdz_fric } else { dpdz_fric };
            let dpdz_psift = dpdz_total_pam / (PSI_TO_PA / FT_TO_M);

            p_est = p_psia + dpdz_psift * seg_len_ft;
        }

        p_psia = p_est;
    }

    p_psia
}

// ============================================================================
//  Gray Gas
// ============================================================================
pub fn gray_fbhp_gas(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> f64 {
    let (tc, pc) = sutton_tc_pc(gsg);
    let osg = 141.5 / (api + 131.5);
    let total_mass = 0.0765 * gsg * qg_mmscfd * 1e6
        + osg * 62.4 * cgr * qg_mmscfd * 5.615
        + wsg * 62.4 * qw_bwpd * 5.615;
    if total_mass < 1e-6 || qg_mmscfd < 0.05 {
        return static_gas_column_pressure(thp, length, tht, bht, gsg, theta);
    }

    let diam_ft = tid / 12.0;
    let rough_ft = rough / 12.0;
    let area = std::f64::consts::PI * diam_ft * diam_ft / 4.0;

    let ndiv = calc_segments(length);
    let seg_len = length / ndiv as f64;

    let m_flow_g_lbs = 0.0765 * gsg * qg_mmscfd * 1e6 / 86400.0;
    let m_flow_w_lbs = wsg * 62.4 * qw_bwpd * 5.615 / 86400.0;

    let mut p_psia = thp;

    for i in 1..=ndiv {
        let frac = (i as f64 - 0.5) / ndiv as f64;
        let temp_f_i = tht + (bht - tht) * frac;
        let mut p_est = p_psia;

        for _ in 0..2 {
            let mut p_avg = (p_psia + p_est) / 2.0;
            if p_avg < 14.7 {
                p_avg = 14.7;
            }

            let (cgr_loc, qo_loc, ql_loc, lsg_loc) =
                condensate_dropout(cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg);

            let m_flow_o_lbs = osg * 62.4 * qo_loc * 5.615 / 86400.0;
            let m_flow_l_lbs = m_flow_o_lbs + m_flow_w_lbs;
            let m_flow_total_lbs = m_flow_g_lbs + m_flow_l_lbs;

            let oil_vis_loc =
                condensate_vis(pr, cgr_loc, gsg, api, temp_f_i, p_avg, oil_vis);

            let zee = z_factor(gsg, temp_f_i, p_avg, tc, pc);
            let mu_g_cp = gas_viscosity(gsg, temp_f_i, p_avg, zee, tc, pc);

            let temp_r = temp_f_i + 459.67;
            let rho_g = MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r);
            let rho_l = lsg_loc * 62.4;

            let v_sg = (m_flow_g_lbs / rho_g.max(1e-10)) / area;
            let v_sl = (m_flow_l_lbs / rho_l.max(1e-10)) / area;
            let v_m = v_sg + v_sl;
            let lambda_l = if v_m > 1e-10 { v_sl / v_m } else { 0.0 };
            let rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l);

            let rs_est = standing_rs(gsg, p_avg, temp_f_i, api);
            let sigma = interfacial_tension(
                p_avg, temp_f_i, api, rs_est,
                m_flow_o_lbs * 0.453592,
                m_flow_w_lbs * 0.453592,
                m_flow_g_lbs * 0.453592,
            );

            let hl = gray_liquid_holdup(v_m, rho_l, rho_g, sigma, diam_ft, lambda_l, p_avg);
            let rho_s = rho_l * hl + rho_g * (1.0 - hl);

            let eps_eff = gray_effective_roughness(rough_ft, sigma, rho_ns, v_m, lambda_l);
            let eps_d = eps_eff / diam_ft;

            let water_visc = water_viscosity(p_avg, temp_f_i, 0.0);
            let mu_l_cp = if ql_loc > 0.0 {
                (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc
            } else {
                oil_vis_loc
            };
            let mu_ns = mu_l_cp * lambda_l + mu_g_cp * (1.0 - lambda_l);
            let mu_ns_lbfts = mu_ns * 6.7197e-4;
            let n_re = if mu_ns_lbfts > 0.0 {
                rho_ns * v_m * diam_ft / mu_ns_lbfts
            } else {
                0.0
            };

            let f_fanning = serghides_fanning(n_re, eps_d);
            let f_moody = 4.0 * f_fanning;

            let gm = m_flow_total_lbs / area;

            let dpdz_hydro = rho_s * theta.sin() / 144.0;
            let dpdz_fric = if rho_ns > 0.0 {
                f_moody * gm * gm / (2.0 * GC * diam_ft * rho_ns * 144.0)
            } else {
                0.0
            };
            let ek = gm * v_sg / (GC * p_avg * 144.0);
            let denom_val = (1.0 - ek).max(0.1);

            let dpdz_total =
                (dpdz_hydro + if injection { -dpdz_fric } else { dpdz_fric }) / denom_val;
            p_est = p_psia + dpdz_total * seg_len;
        }

        p_psia = p_est;
    }

    p_psia
}

// ============================================================================
//  Beggs & Brill Gas
// ============================================================================
pub fn bb_fbhp_gas(
    thp: f64, api: f64, gsg: f64, tid: f64, rough: f64,
    length: f64, tht: f64, bht: f64, wsg: f64,
    qg_mmscfd: f64, cgr: f64, qw_bwpd: f64, oil_vis: f64,
    injection: bool, pr: f64, theta: f64,
) -> f64 {
    let (tc, pc) = sutton_tc_pc(gsg);
    let osg = 141.5 / (api + 131.5);
    let total_mass = 0.0765 * gsg * qg_mmscfd * 1e6
        + osg * 62.4 * cgr * qg_mmscfd * 5.615
        + wsg * 62.4 * qw_bwpd * 5.615;
    if total_mass < 1e-6 || qg_mmscfd < 0.05 {
        return static_gas_column_pressure(thp, length, tht, bht, gsg, theta);
    }

    let diam_ft = tid / 12.0;
    let rough_ft = rough / 12.0;
    let area = std::f64::consts::PI * diam_ft * diam_ft / 4.0;
    let eps_d = rough_ft / diam_ft;

    let ndiv = calc_segments(length);
    let seg_len = length / ndiv as f64;

    let m_flow_g_lbs = 0.0765 * gsg * qg_mmscfd * 1e6 / 86400.0;
    let m_flow_w_lbs = wsg * 62.4 * qw_bwpd * 5.615 / 86400.0;

    let mut p_psia = thp;

    for i in 1..=ndiv {
        let frac = (i as f64 - 0.5) / ndiv as f64;
        let temp_f_i = tht + (bht - tht) * frac;
        let mut p_est = p_psia;

        for _ in 0..2 {
            let mut p_avg = (p_psia + p_est) / 2.0;
            if p_avg < 14.7 {
                p_avg = 14.7;
            }

            let (cgr_loc, qo_loc, ql_loc, lsg_loc) =
                condensate_dropout(cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg);

            let m_flow_o_lbs = osg * 62.4 * qo_loc * 5.615 / 86400.0;
            let m_flow_l_lbs = m_flow_o_lbs + m_flow_w_lbs;

            let oil_vis_loc =
                condensate_vis(pr, cgr_loc, gsg, api, temp_f_i, p_avg, oil_vis);

            let zee = z_factor(gsg, temp_f_i, p_avg, tc, pc);
            let mu_g_cp = gas_viscosity(gsg, temp_f_i, p_avg, zee, tc, pc);

            let temp_r = temp_f_i + 459.67;
            let rho_g = MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r);
            let rho_l = lsg_loc * 62.4;

            let v_sg = (m_flow_g_lbs / rho_g.max(1e-10)) / area;
            let v_sl = (m_flow_l_lbs / rho_l.max(1e-10)) / area;
            let v_m = v_sg + v_sl;
            let lambda_l = if v_m > 1e-10 { v_sl / v_m } else { 0.0 };
            let rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l);

            let froude = v_m * v_m / (G_FT * diam_ft);
            let (pattern, trans_a) = bb_flow_pattern(froude, lambda_l);

            let mut hl0 = if pattern == BB_TRANSITION {
                let hl0_seg = bb_horizontal_holdup(lambda_l, froude, BB_SEGREGATED);
                let hl0_int = bb_horizontal_holdup(lambda_l, froude, BB_INTERMITTENT);
                trans_a * hl0_seg + (1.0 - trans_a) * hl0_int
            } else {
                bb_horizontal_holdup(lambda_l, froude, pattern)
            };

            if pattern == BB_SEGREGATED || pattern == BB_INTERMITTENT || pattern == BB_TRANSITION {
                hl0 *= 0.924;
                hl0 = hl0.max(lambda_l);
            }

            let rs_est = standing_rs(gsg, p_avg, temp_f_i, api);
            let sigma = interfacial_tension(
                p_avg, temp_f_i, api, rs_est,
                m_flow_o_lbs * 0.453592,
                m_flow_w_lbs * 0.453592,
                m_flow_g_lbs * 0.453592,
            );
            let sigma_lbf_ft = sigma * 6.852e-5;
            let n_lv = if sigma_lbf_ft > 0.0 {
                1.938 * v_sl * (rho_l / sigma_lbf_ft).powf(0.25)
            } else {
                0.0
            };

            let hl_theta = if pattern == BB_TRANSITION {
                let hl_seg = bb_inclination_correction(
                    bb_horizontal_holdup(lambda_l, froude, BB_SEGREGATED) * 0.924,
                    lambda_l, n_lv, froude, BB_SEGREGATED, theta,
                );
                let hl_int = bb_inclination_correction(
                    bb_horizontal_holdup(lambda_l, froude, BB_INTERMITTENT) * 0.924,
                    lambda_l, n_lv, froude, BB_INTERMITTENT, theta,
                );
                trans_a * hl_seg + (1.0 - trans_a) * hl_int
            } else {
                bb_inclination_correction(hl0, lambda_l, n_lv, froude, pattern, theta)
            };

            let hl_theta = clamp(hl_theta, lambda_l, 1.0);
            let rho_s = rho_l * hl_theta + rho_g * (1.0 - hl_theta);

            let water_visc = water_viscosity(p_avg, temp_f_i, 0.0);
            let mu_l_cp = if ql_loc > 0.0 {
                (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc
            } else {
                oil_vis_loc
            };
            let mu_ns = mu_l_cp * lambda_l + mu_g_cp * (1.0 - lambda_l);
            let mu_ns_lbfts = mu_ns * 6.7197e-4;
            let n_re = if mu_ns_lbfts > 0.0 {
                rho_ns * v_m * diam_ft / mu_ns_lbfts
            } else {
                0.0
            };

            let f_ns = serghides_fanning(n_re, eps_d);
            let f_tp = bb_two_phase_friction(f_ns, lambda_l, hl_theta);

            let dpdz_hydro = rho_s * theta.sin() / 144.0;
            let dpdz_fric = if rho_ns > 0.0 {
                4.0 * f_tp * rho_ns * v_m * v_m / (2.0 * GC * diam_ft * 144.0)
            } else {
                0.0
            };
            let ek = rho_s * v_m * v_sg / (GC * p_avg * 144.0);
            let denom_val = (1.0 - ek).max(0.1);

            let dpdz_total =
                (dpdz_hydro + if injection { -dpdz_fric } else { dpdz_fric }) / denom_val;
            p_est = p_psia + dpdz_total * seg_len;
        }

        p_psia = p_est;
    }

    p_psia
}
