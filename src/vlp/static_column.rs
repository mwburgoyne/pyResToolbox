/// Static column pressure helpers.
/// Direct port of nodal.py lines 685-721.

use super::pvt_helpers::*;

const MW_AIR: f64 = 28.97;

/// Static gas column pressure (psia).
pub fn static_gas_column_pressure(
    thp: f64, length: f64, tht: f64, bht: f64, gsg: f64, theta: f64,
) -> f64 {
    let (tc, pc) = sutton_tc_pc(gsg);
    let n_seg = 50;
    let d_len = length / n_seg as f64;
    let sin_theta = theta.sin();
    let mut p = thp;
    for i in 0..n_seg {
        let frac = (i as f64 + 0.5) / n_seg as f64;
        let temp_local = tht + (bht - tht) * frac;
        let temp_r = temp_local + 459.67;
        let zee = z_factor(gsg, temp_local, p.max(14.7), tc, pc);
        let rho_gas = MW_AIR * gsg * p / (zee * 10.732 * temp_r);
        p += rho_gas / 144.0 * d_len * sin_theta;
    }
    p
}

/// Static oil column pressure (psia).
pub fn static_oil_column_pressure(
    thp: f64, length: f64, tht: f64, bht: f64, wc: f64, wsg: f64,
    api: f64, sgsp: f64, pb: f64, rsb: f64, rsb_scale: f64, theta: f64,
) -> f64 {
    let sgsto = 141.5 / (api + 131.5);
    let rsb_for_calc = rsb / rsb_scale;
    let sin_theta = theta.sin();
    let n_seg = 50;
    let d_len = length / n_seg as f64;
    let mut p = thp;
    for i in 0..n_seg {
        let frac = (i as f64 + 0.5) / n_seg as f64;
        let temp_local = tht + (bht - tht) * frac;
        let rho_oil = if p < pb {
            let rs = velarde_rs(sgsp, api, temp_local, pb, rsb_for_calc, p) * rsb_scale;
            oil_density_mccain(rs, sgsp, sgsto, p, temp_local)
        } else {
            oil_density_mccain(rsb, sgsp, sgsto, pb, temp_local)
        };
        let oil_sg_local = rho_oil / 62.4;
        let mix_sg = (1.0 - wc) * oil_sg_local + wc * wsg;
        p += 0.433 * mix_sg * d_len * sin_theta;
    }
    p
}
