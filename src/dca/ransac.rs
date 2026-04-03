/// RANSAC-based linear regression with OLS fallback.
///
/// Parameters:
///   x, y        — data arrays (same length)
///   n_iter       — number of RANSAC iterations
///   threshold_sigma — inlier threshold in multiples of MAD-estimated sigma
///   seed         — seed for deterministic PRNG
///   through_origin — if true, fit y = slope*x (intercept forced to 0)
///
/// Returns (slope, intercept, inlier_mask).

/// Simple Park-Miller LCG PRNG (state must be > 0).
struct ParkMillerRng {
    state: u64,
}

impl ParkMillerRng {
    fn new(seed: u64) -> Self {
        // Ensure state is never zero
        let s = if seed == 0 { 1 } else { seed };
        ParkMillerRng { state: s }
    }

    fn next_u64(&mut self) -> u64 {
        // Park-Miller: next = (state * 48271) mod 2147483647
        self.state = (self.state.wrapping_mul(48271)) % 2_147_483_647;
        self.state
    }

    /// Return a random index in [0, n).
    fn rand_index(&mut self, n: usize) -> usize {
        (self.next_u64() % (n as u64)) as usize
    }

    /// Choose 2 distinct indices from [0, n) without replacement.
    fn choice2(&mut self, n: usize) -> (usize, usize) {
        let i1 = self.rand_index(n);
        let mut i2 = self.rand_index(n - 1);
        if i2 >= i1 {
            i2 += 1;
        }
        (i1, i2)
    }
}

/// Compute median of a slice (sorts a copy).
fn median(vals: &[f64]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    let mut sorted = vals.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
    }
}

/// Ordinary least-squares fit.
/// If through_origin: y = slope * x  (intercept = 0).
/// Otherwise: y = slope * x + intercept.
fn ols(x: &[f64], y: &[f64], through_origin: bool) -> (f64, f64) {
    let n = x.len() as f64;
    if through_origin {
        let mut sum_xy = 0.0;
        let mut sum_xx = 0.0;
        for i in 0..x.len() {
            sum_xy += x[i] * y[i];
            sum_xx += x[i] * x[i];
        }
        let slope = if sum_xx.abs() < 1e-30 { 0.0 } else { sum_xy / sum_xx };
        (slope, 0.0)
    } else {
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut sum_xy = 0.0;
        let mut sum_xx = 0.0;
        for i in 0..x.len() {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            sum_xx += x[i] * x[i];
        }
        let denom = n * sum_xx - sum_x * sum_x;
        if denom.abs() < 1e-30 {
            let intercept = sum_y / n;
            (0.0, intercept)
        } else {
            let slope = (n * sum_xy - sum_x * sum_y) / denom;
            let intercept = (sum_y - slope * sum_x) / n;
            (slope, intercept)
        }
    }
}

/// OLS fit on a subset of points indicated by a boolean mask.
fn ols_masked(x: &[f64], y: &[f64], mask: &[bool], through_origin: bool) -> (f64, f64) {
    let xsub: Vec<f64> = x.iter().zip(mask).filter(|(_, &m)| m).map(|(&v, _)| v).collect();
    let ysub: Vec<f64> = y.iter().zip(mask).filter(|(_, &m)| m).map(|(&v, _)| v).collect();
    ols(&xsub, &ysub, through_origin)
}

/// OLS fit on exactly two points given by indices.
fn ols_two(x: &[f64], y: &[f64], i1: usize, i2: usize, through_origin: bool) -> (f64, f64) {
    if through_origin {
        let sum_xy = x[i1] * y[i1] + x[i2] * y[i2];
        let sum_xx = x[i1] * x[i1] + x[i2] * x[i2];
        let slope = if sum_xx.abs() < 1e-30 { 0.0 } else { sum_xy / sum_xx };
        (slope, 0.0)
    } else {
        let dx = x[i2] - x[i1];
        if dx.abs() < 1e-30 {
            let intercept = 0.5 * (y[i1] + y[i2]);
            (0.0, intercept)
        } else {
            let slope = (y[i2] - y[i1]) / dx;
            let intercept = y[i1] - slope * x[i1];
            (slope, intercept)
        }
    }
}

pub fn ransac_linreg(
    x: &[f64],
    y: &[f64],
    n_iter: usize,
    threshold_sigma: f64,
    seed: u64,
    through_origin: bool,
) -> (f64, f64, Vec<bool>) {
    let n = x.len();
    assert_eq!(n, y.len(), "x and y must have the same length");

    // Full-data OLS as fallback
    let (full_slope, full_intercept) = ols(x, y, through_origin);
    let all_inliers = vec![true; n];

    if n < 3 {
        return (full_slope, full_intercept, all_inliers);
    }

    // Compute residuals from full OLS
    let residuals: Vec<f64> = (0..n)
        .map(|i| y[i] - (full_slope * x[i] + full_intercept))
        .collect();

    // MAD-based sigma estimate
    let med_res = median(&residuals);
    let abs_devs: Vec<f64> = residuals.iter().map(|&r| (r - med_res).abs()).collect();
    let mad = median(&abs_devs);
    let sigma_est = 1.4826 * mad;

    // If sigma is effectively zero, all points are perfectly on the line
    if sigma_est < 1e-15 {
        return (full_slope, full_intercept, all_inliers);
    }

    let threshold = threshold_sigma * sigma_est;
    let min_samples = 2_usize;
    let min_inliers = min_samples.max(((0.5 * n as f64).ceil()) as usize);

    let mut rng = ParkMillerRng::new(seed);
    let mut best_n_inliers = 0_usize;
    let mut best_slope = full_slope;
    let mut best_intercept = full_intercept;
    let mut best_mask = all_inliers.clone();

    for _ in 0..n_iter {
        // Pick 2 random distinct points
        let (i1, i2) = rng.choice2(n);

        // Fit on subset
        let (slope, intercept) = ols_two(x, y, i1, i2, through_origin);

        // Count inliers
        let mut inlier_mask = vec![false; n];
        let mut n_inliers = 0_usize;
        for i in 0..n {
            let res = (y[i] - (slope * x[i] + intercept)).abs();
            if res < threshold {
                inlier_mask[i] = true;
                n_inliers += 1;
            }
        }

        if n_inliers > best_n_inliers && n_inliers >= min_inliers {
            // Refit on all inliers
            let (refit_slope, refit_intercept) = ols_masked(x, y, &inlier_mask, through_origin);
            best_slope = refit_slope;
            best_intercept = refit_intercept;
            best_n_inliers = n_inliers;
            best_mask = inlier_mask;
        }
    }

    // If we never found enough inliers, fall back to full OLS
    if best_n_inliers < min_inliers {
        return (full_slope, full_intercept, vec![true; n]);
    }

    (best_slope, best_intercept, best_mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_line() {
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 2.0 * xi + 3.0).collect();
        let (slope, intercept, mask) = ransac_linreg(&x, &y, 100, 3.0, 42, false);
        assert!((slope - 2.0).abs() < 1e-10);
        assert!((intercept - 3.0).abs() < 1e-10);
        assert!(mask.iter().all(|&m| m));
    }

    #[test]
    fn test_through_origin() {
        let x: Vec<f64> = (1..11).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 5.0 * xi).collect();
        let (slope, intercept, _mask) = ransac_linreg(&x, &y, 100, 3.0, 42, true);
        assert!((slope - 5.0).abs() < 1e-10);
        assert!(intercept.abs() < 1e-10);
    }
}
