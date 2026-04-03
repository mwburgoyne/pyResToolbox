/// Material balance objective: computes CV (std/mean) of OOIP estimates.
/// Direct port of matbal.py lines 592-609.

pub fn matbal_objective(
    params: &[f64],
    param_names: &[String],
    f_arr: &[f64],
    eo_arr: &[f64],
    eg_arr: &[f64],
    p_field: &[f64],
    boi: f64,
    base_m: f64,
    base_cf: f64,
    base_cw: f64,
    base_sw_i: f64,
) -> f64 {
    let mut m = base_m;
    let mut cf = base_cf;
    let mut cw = base_cw;
    let mut sw_i = base_sw_i;

    for (k, v) in param_names.iter().zip(params.iter()) {
        match k.as_str() {
            "m" => m = *v,
            "cf" => cf = *v,
            "cw" => cw = *v,
            "sw_i" => sw_i = *v,
            _ => {}
        }
    }

    let n = f_arr.len();
    let mut efw = vec![0.0f64; n];
    for i in 0..n {
        let dp = p_field[0] - p_field[i];
        if (1.0 - sw_i) > 0.0 {
            efw[i] = boi * (cw * sw_i + cf) / (1.0 - sw_i) * dp;
        }
    }

    let mut n_estimates: Vec<f64> = Vec::new();
    for i in 0..n {
        let denom = eo_arr[i] + m * eg_arr[i] + (1.0 + m) * efw[i];
        if denom.abs() > 1e-30 {
            n_estimates.push(f_arr[i] / denom);
        }
    }

    if n_estimates.is_empty() {
        return 1e10;
    }

    let mean_n: f64 = n_estimates.iter().sum::<f64>() / n_estimates.len() as f64;
    if mean_n.abs() < 1e-30 {
        return 1e10;
    }

    let variance: f64 = n_estimates.iter().map(|x| (x - mean_n).powi(2)).sum::<f64>() / n_estimates.len() as f64;
    let std_dev = variance.sqrt();
    std_dev / mean_n.abs()
}
