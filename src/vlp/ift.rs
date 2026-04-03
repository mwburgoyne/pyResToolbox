/// Interfacial tension correlations for VLP.
/// Direct port of nodal.py lines 589-628.

use super::pvt_helpers::clamp;

const SCF_STB_TO_M3M3: f64 = 0.178108;

/// Baker & Swerdloff dead oil IFT (dyne/cm).
fn dead_oil_ift(api: f64, temp_f: f64) -> f64 {
    let temp_c = (temp_f - 32.0) / 1.8;
    let a = 1.11591 - 0.00305 * temp_c;
    (a * (38.085 - 0.259 * api)).max(1.0)
}

/// Firoozabadi & Ramey gas-oil IFT (dyne/cm).
pub fn gas_oil_ift(api: f64, temp_f: f64, rs_scf_stb: f64) -> f64 {
    let sigma_od = dead_oil_ift(api, temp_f);
    if rs_scf_stb <= 0.0 {
        return sigma_od;
    }
    let rs_m3m3 = rs_scf_stb * SCF_STB_TO_M3M3;
    let ratio = if rs_m3m3 < 50.0 {
        1.0 / (1.0 + 0.02549 * rs_m3m3.powf(1.0157))
    } else {
        32.0436 * rs_m3m3.powf(-1.1367)
    };
    (sigma_od * clamp(ratio, 0.0, 1.0)).max(1.0)
}

/// Jennings & Newman gas-water IFT (dyne/cm).
pub fn gas_water_ift(press_psia: f64, temp_f: f64) -> f64 {
    let p = press_psia.max(14.7);
    let sw74 = 75.0 - 1.108 * p.powf(0.349);
    let sw280 = 53.0 - 0.1048 * p.powf(0.637);
    if temp_f <= 74.0 {
        return sw74.max(1.0);
    } else if temp_f >= 280.0 {
        return sw280.max(1.0);
    }
    (sw74 + (sw280 - sw74) * (temp_f - 74.0) / (280.0 - 74.0)).max(1.0)
}

/// Blended IFT based on mass flow fractions.
pub fn interfacial_tension(
    p_avg: f64, temp_f: f64, api: f64, rs_scf_stb: f64,
    m_flow_o: f64, m_flow_w: f64, m_flow_g: f64,
) -> f64 {
    let sigma_og = gas_oil_ift(api, temp_f, rs_scf_stb);
    let sigma_wg = gas_water_ift(p_avg, temp_f);
    let sigma_ow = 26.0;
    let denom = m_flow_o * m_flow_w + m_flow_o * m_flow_g + m_flow_g * m_flow_w;
    if denom <= 0.0 {
        return 20.0;
    }
    (((m_flow_o * m_flow_w) * sigma_ow
        + (m_flow_o * m_flow_g) * sigma_og
        + (m_flow_g * m_flow_w) * sigma_wg) / denom)
    .max(1.0)
}
