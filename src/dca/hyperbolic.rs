use crate::dca::ransac::ransac_linreg;

/// Fit hyperbolic decline q = qi / (1 + b*di*t)^(1/b) to rate-vs-time data.
///
/// Grid search over b in [0.05, 0.95] step 0.01.
/// For each b, linearize: Y = q^(-b) = qi^(-b) + qi^(-b)*b*di * t
/// Fit Y vs t with RANSAC, extract qi and di, compute R2 in original rate space.
///
/// Returns Some((qi, di, b, R2)) or None if no valid fit found.
pub fn fit_hyperbolic(t: &[f64], q: &[f64]) -> Option<(f64, f64, f64, f64)> {
    let n = t.len();
    if n < 3 {
        return None;
    }

    // Precompute ss_tot for R2
    let q_mean: f64 = q.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = q.iter().map(|&qi| (qi - q_mean).powi(2)).sum();
    if ss_tot < 1e-30 {
        return None;
    }

    let mut best_r2 = f64::NEG_INFINITY;
    let mut best_result: Option<(f64, f64, f64, f64)> = None;

    // Grid search: b from 0.05 to 0.95 in steps of 0.01 (91 values)
    for ib in 5..=95 {
        let b = ib as f64 * 0.01;

        // Linearize: Y = q^(-b)
        let mut y_lin: Vec<f64> = Vec::with_capacity(n);
        let mut valid = true;
        for i in 0..n {
            if q[i] <= 0.0 {
                valid = false;
                break;
            }
            y_lin.push(q[i].powf(-b));
        }
        if !valid {
            continue;
        }

        // RANSAC fit: Y = intercept + slope * t
        let (slope, intercept, _mask) = ransac_linreg(t, &y_lin, 200, 3.0, 42, false);

        // intercept = qi^(-b), must be positive
        if intercept <= 0.0 {
            continue;
        }

        // qi = intercept^(-1/b)
        let qi = intercept.powf(-1.0 / b);
        if qi <= 0.0 || !qi.is_finite() {
            continue;
        }

        // slope = qi^(-b) * b * di  =>  di = slope / (intercept * b)
        let di = slope / (intercept * b);
        if di <= 0.0 || !di.is_finite() {
            continue;
        }

        // Compute predicted rates and R2
        let inv_b = 1.0 / b;
        let mut ss_res = 0.0;
        for i in 0..n {
            let inner = 1.0 + b * di * t[i];
            let q_pred = if inner > 0.0 {
                qi / inner.powf(inv_b)
            } else {
                0.0
            };
            ss_res += (q[i] - q_pred).powi(2);
        }

        let r2 = 1.0 - ss_res / ss_tot;

        if r2 > best_r2 {
            best_r2 = r2;
            best_result = Some((qi, di, b, r2));
        }
    }

    best_result
}

/// Fit hyperbolic decline to cumulative production vs rate data.
///
/// Relationship: Np = qi/((1-b)*di) - qi^b / ((1-b)*di) * q^(1-b)
/// Linearize:    Np = A + B * q^(1-b)
///   where A = qi/((1-b)*di), B = -qi^b / ((1-b)*di)
///   ratio B/A = -1/qi^(1-b)  =>  qi = (-A/B)^(1/(1-b))
///   di = qi / ((1-b) * A)
///
/// Returns Some((qi, di, b, R2)) or None if no valid fit found.
pub fn fit_hyperbolic_cum(np_cum: &[f64], q: &[f64]) -> Option<(f64, f64, f64, f64)> {
    let n = np_cum.len();
    if n < 3 {
        return None;
    }

    // Precompute ss_tot for R2 in rate space
    let q_mean: f64 = q.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = q.iter().map(|&qi| (qi - q_mean).powi(2)).sum();
    if ss_tot < 1e-30 {
        return None;
    }

    let mut best_r2 = f64::NEG_INFINITY;
    let mut best_result: Option<(f64, f64, f64, f64)> = None;

    // Grid search: b from 0.05 to 0.95 in steps of 0.01 (91 values)
    for ib in 5..=95 {
        let b = ib as f64 * 0.01;
        let exp = 1.0 - b; // exponent for q

        // Linearize: X = q^(1-b), Y = Np
        // Np = intercept + slope * X
        let mut x_lin: Vec<f64> = Vec::with_capacity(n);
        let mut valid = true;
        for i in 0..n {
            if q[i] <= 0.0 {
                valid = false;
                break;
            }
            x_lin.push(q[i].powf(exp));
        }
        if !valid {
            continue;
        }

        // RANSAC fit: Np = intercept + slope * q^(1-b)
        let (slope, intercept, _mask) = ransac_linreg(&x_lin, np_cum, 200, 3.0, 42, false);

        // ratio = slope / intercept = -1 / qi^(1-b)
        if intercept.abs() < 1e-30 {
            continue;
        }
        let ratio = slope / intercept;

        // ratio must be negative for valid solution
        if ratio >= 0.0 {
            continue;
        }

        // qi^(1-b) = -1/ratio = -intercept/slope
        let qi_exp = -1.0 / ratio; // = -intercept / slope, should be positive
        if qi_exp <= 0.0 {
            continue;
        }

        let qi = qi_exp.powf(1.0 / exp);
        if qi <= 0.0 || !qi.is_finite() {
            continue;
        }

        // di = qi / (exp * intercept)
        // intercept = A = qi / ((1-b)*di), so di = qi / ((1-b)*A)
        let di = qi / (exp * intercept);
        if di <= 0.0 || !di.is_finite() {
            continue;
        }

        // Compute predicted rates from Np:
        // q = qi * (1 - (1-b)*di*Np/qi)^(1/(1-b))
        // equivalently: inner = 1 - exp*di*Np/qi
        //               q_pred = qi * inner^(1/exp)
        let inv_exp = 1.0 / exp;
        let mut ss_res = 0.0;
        for i in 0..n {
            let inner = 1.0 - exp * di * np_cum[i] / qi;
            let inner_clamped = if inner < 1e-10 { 1e-10 } else { inner };
            let q_pred = qi * inner_clamped.powf(inv_exp);
            ss_res += (q[i] - q_pred).powi(2);
        }

        let r2 = 1.0 - ss_res / ss_tot;

        if r2 > best_r2 {
            best_r2 = r2;
            best_result = Some((qi, di, b, r2));
        }
    }

    best_result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_hyperbolic_known() {
        // Generate synthetic data: qi=1000, di=0.05, b=0.5
        let qi_true = 1000.0;
        let di_true = 0.05;
        let b_true = 0.5;
        let inv_b = 1.0 / b_true;

        let t: Vec<f64> = (0..60).map(|i| (i + 1) as f64).collect();
        let q: Vec<f64> = t
            .iter()
            .map(|&ti| qi_true / (1.0 + b_true * di_true * ti).powf(inv_b))
            .collect();

        let result = fit_hyperbolic(&t, &q);
        assert!(result.is_some());
        let (qi, di, b, r2) = result.unwrap();
        assert!((qi - qi_true).abs() / qi_true < 0.02, "qi={}", qi);
        assert!((di - di_true).abs() / di_true < 0.02, "di={}", di);
        assert!((b - b_true).abs() < 0.02, "b={}", b);
        assert!(r2 > 0.999, "r2={}", r2);
    }

    #[test]
    fn test_fit_hyperbolic_cum_known() {
        // Generate synthetic data: qi=1000, di=0.05, b=0.5
        let qi_true = 1000.0;
        let di_true = 0.05;
        let b_true = 0.5;
        let inv_b = 1.0 / b_true;
        let exp = 1.0 - b_true;

        let t: Vec<f64> = (0..60).map(|i| (i + 1) as f64).collect();
        let q: Vec<f64> = t
            .iter()
            .map(|&ti| qi_true / (1.0 + b_true * di_true * ti).powf(inv_b))
            .collect();

        // Compute cumulative production analytically:
        // Np(t) = qi / ((1-b)*di) * (1 - (1+b*di*t)^(-(1-b)/b))
        //       = qi / ((1-b)*di) * (1 - (1+b*di*t)^(-exp/b))
        // Actually: Np = qi/((1-b)*di) * [1 - (1+b*di*t)^((1-1/b))]
        // Np = qi/((1-b)*di) * [1 - 1/(1+b*di*t)^((1-b)/b)]
        let np_cum: Vec<f64> = t
            .iter()
            .map(|&ti| {
                let factor = (1.0 + b_true * di_true * ti).powf((1.0 - b_true) / b_true);
                qi_true / (exp * di_true) * (1.0 - 1.0 / factor)
            })
            .collect();

        let result = fit_hyperbolic_cum(&np_cum, &q);
        assert!(result.is_some());
        let (qi, di, b, r2) = result.unwrap();
        assert!((qi - qi_true).abs() / qi_true < 0.02, "qi={}", qi);
        assert!((di - di_true).abs() / di_true < 0.02, "di={}", di);
        assert!((b - b_true).abs() < 0.02, "b={}", b);
        assert!(r2 > 0.999, "r2={}", r2);
    }
}
