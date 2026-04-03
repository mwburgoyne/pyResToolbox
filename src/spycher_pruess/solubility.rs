// =========================================================================
// Spycher-Pruess CO2-Brine mutual solubility model
//
// Port of CO2_Brine_Mixture.co2BrineSolubility() from pyrestoolbox/brine/brine.py
//
// Reference:
//   Spycher & Pruess, "A Phase-Partitioning Model for CO2-Brine Mixtures
//   at Elevated Temperatures and Pressures: Application to CO2-Enhanced
//   Geothermal Systems", Transp Porous Med (2010) 82:173-196
// =========================================================================

// ---- Physical constants ------------------------------------------------
const MWSAL: f64 = 58.4428;       // Mole weight of NaCl
const MWWAT: f64 = 18.01528;      // Mole weight of pure water
const MWCO2: f64 = 44.01;         // Mole weight of CO2
const RGASCON: f64 = 83.1447;     // bar.cm3/(mol.K)
const CEL2KEL: f64 = 273.15;      // degK at 0 degC
const CONMOLA: f64 = 1000.0 / MWWAT; // Moles in 1000 kg water (molality conversion)
const BAR2PSI: f64 = 14.5037738;
const EPS: f64 = 1e-8;

// ---- Helper: polynomial evaluator  sum(x[i] * t^i) --------------------
#[inline]
fn ft(t: f64, coeffs: &[f64]) -> f64 {
    let mut result = 0.0;
    let mut t_pow = 1.0;
    for &c in coeffs {
        result += c * t_pow;
        t_pow *= t;
    }
    result
}

// ---- State carried through the calculation -----------------------------
struct SpState {
    p_bar: f64,
    deg_c: f64,
    t_kel: f64,
    ppm: f64,
    p_rt: f64,   // pBar / (R * T)
    p_rt0: f64,  // (pBar - P0) / (R * T)
    p0: f64,

    molal: f64,   // brine molality (gmol/kg)
    x_salt: f64,  // xNa + xCl in aqueous phase

    y: [f64; 2],  // [yCO2, yH2O] in gas phase
    x: [f64; 2],  // [xCO2, xH2O] in aqueous phase

    gamma_prime: f64,
    gamma: [f64; 2],

    // RK-EOS parameters
    a_mix: f64,
    aij: [[f64; 2]; 2],
    kij: [[f64; 2]; 2],
    b_mix: f64,
    b: [f64; 2],

    molar_vol: f64, // cm3/gmol
    v_bar: [f64; 2], // avg partial molar volume of condensed phase

    fug_pi: [f64; 2], // fugacity * pressure
    k_val: [f64; 2],  // K-values (yi/xi)

    a_eq: f64,       // equilibrium factor A (Eq 10)
    b_prime: f64,    // equilibrium factor B' (Eq 17)

    co2_sat: bool,   // whether CO2 is saturated liquid
    low_temp: bool,  // calculation path flag
    scaled: bool,    // blended 99-109 degC path flag
}

impl SpState {
    fn new(p_bar: f64, deg_c: f64, ppm: f64) -> Self {
        let t_kel = deg_c + CEL2KEL;
        SpState {
            p_bar,
            deg_c,
            t_kel,
            ppm,
            p_rt: 0.0,
            p_rt0: 0.0,
            p0: 0.0,
            molal: 0.0,
            x_salt: 0.0,
            y: [1.0, 0.0],
            x: [0.0, 1.0],
            gamma_prime: 1.0,
            gamma: [1.0, 1.0],
            a_mix: 0.0,
            aij: [[0.0; 2]; 2],
            kij: [[0.0; 2]; 2],
            b_mix: 0.0,
            b: [0.0; 2],
            molar_vol: 0.0,
            v_bar: [0.0; 2],
            fug_pi: [0.0; 2],
            k_val: [0.0; 2],
            a_eq: 0.0,
            b_prime: 0.0,
            co2_sat: false,
            low_temp: true,
            scaled: false,
        }
    }

    // ---- ppm -> molality (gmol NaCl / kg water) -----------------------
    fn ppm2molality(&self) -> f64 {
        self.ppm / MWSAL * 1000.0 / (1e6 - self.ppm)
    }

    // ---- Blended value between 99-109 degC ----------------------------
    #[inline]
    fn blended_val(&self, low_val: f64, high_val: f64) -> f64 {
        ((self.deg_c - 99.0) * low_val + (109.0 - self.deg_c) * high_val) / 10.0
    }

    // ---- Determine calculation type -----------------------------------
    fn calc_type(&mut self) {
        self.low_temp = true;
        self.scaled = false;
        if self.deg_c > 99.0 && self.deg_c <= 109.0 {
            self.low_temp = false;
            self.scaled = true;
        } else if self.deg_c > 109.0 {
            self.low_temp = false;
            self.scaled = false;
        }
    }

    // ---- Water vapor pressure (Buck equation) -------------------------
    fn water_vap_p(&self) -> f64 {
        let kpa = 0.61121
            * ((18.678 - self.deg_c / 234.5) * (self.deg_c / (257.14 + self.deg_c))).exp();
        0.145038 * kpa // psia
    }

    // ---- Initial estimate for yH2O -----------------------------------
    fn est_yh2o(&self) -> f64 {
        let deg_f = self.deg_c * 1.8 + 32.0;
        let p = self.p_bar * BAR2PSI;

        // Slope: Bleasdale (Shifted power) + c
        let slope = {
            let (a, b, c, d) = (1.939e0, 8.913e-3, -4.844e0, 2.551e-4);
            (a + b * deg_f).powf(c) + d
        };

        // Intercept: Inv pow
        let intcpt = {
            let (a, b, c) = (3.097e-1, 9.136e-6, 1.719e0);
            -1.0 / (a + b * deg_f * deg_f).powf(c)
        };

        let pp_ratio = slope * p + intcpt + 1.0;
        let atm_pvap = self.water_vap_p();
        let mut pp = pp_ratio * atm_pvap;
        if pp < atm_pvap {
            pp = atm_pvap;
        }
        (pp / p).clamp(0.0, 1.0 - EPS)
    }

    // ---- CO2 activity coefficient (gamma prime, Eq 18) ----------------
    fn gamma_co2(&mut self) {
        let cl = [0.0002217, 1.074, 2648.0];
        let cz = [0.000013, -20.12, 5259.0];

        let t = self.t_kel;
        let lamb = cl[0] * t + cl[1] / t + cl[2] / (t * t);
        let zeta = cz[0] * t + cz[1] / t + cz[2] / (t * t);

        self.gamma_prime = (1.0 + self.molal / CONMOLA)
            * (self.molal * (2.0 * lamb + self.molal * zeta)).exp();
    }

    // ---- aCO2 for RK-EOS (temperature dependent) ---------------------
    fn a_co2_rk(&self) -> f64 {
        if self.low_temp {
            ft(self.t_kel, &[7.54e7, -4.13e4]) // Low temp
        } else {
            ft(self.t_kel, &[8.008e7, -4.984e4]) // High temp
        }
    }

    // ---- aH2O for RK-EOS (only used high temp) -----------------------
    fn a_h2o_rk(&self) -> f64 {
        ft(self.t_kel, &[1.337e8, -1.4e4])
    }

    // ---- aMix for RK-EOS (mixing rules) ------------------------------
    fn a_mix_rk(&mut self) {
        // K array (big K, not kij)
        let mut k_big = [[0.0_f64; 2]; 2];
        if !self.low_temp {
            k_big[0][1] = ft(self.t_kel, &[0.4228, -7.422e-4]); // K_CO2-H2O
            k_big[1][0] = ft(self.t_kel, &[1.427e-2, -4.037e-4]); // K_H2O-CO2
        } else {
            k_big[0][1] = 7.89e7;
            k_big[1][0] = 7.89e7;
        }

        // kij array (Eq A-6)
        let mut kij = [[0.0_f64; 2]; 2];
        kij[0][1] = k_big[0][1] * self.y[0] + k_big[1][0] * self.y[1];
        kij[1][0] = kij[0][1];

        // aij array
        let mut aij = [[0.0_f64; 2]; 2];
        aij[0][0] = self.a_co2_rk();
        aij[1][1] = self.a_h2o_rk();
        if self.low_temp {
            aij[0][1] = 7.89e7;
            aij[1][0] = 7.89e7;
        } else {
            for i in 0..2_usize {
                let j = 1 - i;
                aij[i][j] = (aij[i][i] * aij[j][j]).sqrt() * (1.0 - kij[i][j]); // Eq A-5
            }
        }

        // amix = sum_i sum_j y[i]*y[j]*aij[i][j]
        let mut amix = 0.0;
        for i in 0..2 {
            for j in 0..2 {
                amix += self.y[i] * self.y[j] * aij[i][j];
            }
        }

        self.a_mix = amix;
        self.aij = aij;
        self.kij = kij;
    }

    // ---- bMix for RK-EOS ---------------------------------------------
    fn b_mix_rk(&mut self) {
        let b = if self.low_temp {
            [27.80, 18.18]
        } else {
            [28.25, 15.70]
        };
        self.b_mix = self.y[0] * b[0] + self.y[1] * b[1];
        self.b = b;
    }

    // ---- Cubic solver (Halley's method, all real roots) ---------------
    // Solves V^3 + e2*V^2 + e1*V + e0 = 0, returns molar volume
    // Selects root via Gibbs energy criterion (Eqs 25-26, Spycher & Pruess 2003)
    fn cubic_solver(&mut self, e2: f64, e1: f64, e0: f64) -> f64 {
        let roots = halley_all_roots(e2, e1, e0);

        if roots.len() > 1 {
            let vgas = roots.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let vliq = roots.iter().cloned().fold(f64::INFINITY, f64::min);

            let w1 = self.p_bar * (vgas - vliq);

            // Gibbs energy comparison
            let ln_arg1 = if (vgas - self.b_mix).abs() < 1e-30 {
                1e-30
            } else {
                (vgas - self.b_mix) / (vliq - self.b_mix)
            };
            let ln_arg2 = if ((vliq + self.b_mix) * vgas).abs() < 1e-30 {
                1e-30
            } else {
                (vgas + self.b_mix) * vliq / ((vliq + self.b_mix) * vgas)
            };
            let w2 = RGASCON * self.t_kel * ln_arg1.abs().ln()
                + self.a_mix / (self.t_kel.sqrt() * self.b_mix) * ln_arg2.abs().ln();

            if w2 - w1 > 0.0 {
                // Vapor root preferred
                let old_co2_sat = self.co2_sat;
                self.co2_sat = false;
                if old_co2_sat {
                    // CO2 was saturated, now it's not -- need to repeat
                    // (caller checks co2_sat change via repeat flag)
                }
                vgas
            } else {
                // Liquid root preferred
                let old_co2_sat = self.co2_sat;
                self.co2_sat = true;
                if !old_co2_sat {
                    // CO2 was not saturated, now it is
                }
                vliq
            }
        } else if !roots.is_empty() {
            let old_co2_sat = self.co2_sat;
            if old_co2_sat {
                self.co2_sat = false;
            }
            roots[0]
        } else {
            // Fallback: shouldn't happen, but return estimate
            -e2 / 3.0
        }
    }

    // ---- Molar volume from RK-EOS cubic ------------------------------
    fn molar_volume(&mut self) -> bool {
        // Returns true if co2_sat state changed (repeat needed)
        let old_co2_sat = self.co2_sat;

        let rt_p = RGASCON * self.t_kel / self.p_bar;
        let at12p = self.a_mix / (self.p_bar * self.t_kel.sqrt());

        // Coefficients of RK cubic (Eq A-2)
        let e2 = -rt_p;
        let e1 = -(rt_p * self.b_mix - at12p + self.b_mix * self.b_mix);
        let e0 = -at12p * self.b_mix;

        self.molar_vol = self.cubic_solver(e2, e1, e0);

        old_co2_sat != self.co2_sat
    }

    // ---- Fugacity coefficient * pressure (Eq A-8) --------------------
    fn fug_p(&mut self) {
        let y = self.y;
        let x = self.x;
        let kij = self.kij;
        let a_mix = self.a_mix;
        let aij = self.aij;
        let b_mix = self.b_mix;
        let b = self.b;
        let v_mol = self.molar_vol;

        // Clamp compositions
        let y_co2 = y[0].clamp(0.0, 1.0);
        let x_co2 = x[0].clamp(0.0, 1.0);
        let y_arr = [y_co2, 1.0 - y_co2];
        let x_arr = [x_co2, 1.0 - x_co2];
        self.y = y_arr;
        self.x = x_arr;

        for k in 0..2_usize {
            let t1 = b[k] / b_mix * (self.p_rt * v_mol - 1.0);
            let t2 = -(self.p_rt * (v_mol - b_mix).max(1e-9)).ln();

            // t3 computation (Eq A-8 numerator terms)
            let mut t3 = 0.0;
            // First sum: sum_i y[i] * (aij[i][k] + aij[k][i])
            for i in 0..2 {
                t3 += y_arr[i] * (aij[i][k] + aij[k][i]);
            }
            // Subtract double loop: sum_i sum_j y[i]^2 * y[j] * (kij[i][j] - kij[j][i]) * sqrt(aij[i][i]*aij[j][j])
            for i in 0..2 {
                for j in 0..2 {
                    t3 -= y_arr[i] * y_arr[i] * y_arr[j]
                        * (kij[i][j] - kij[j][i])
                        * (aij[i][i] * aij[j][j]).sqrt();
                }
            }
            // Add: sum_i x[k] * x[i] * (kij[k][i] - kij[i][k]) * sqrt(aij[i][i]*aij[k][k])
            for i in 0..2 {
                t3 += x_arr[k] * x_arr[i]
                    * (kij[k][i] - kij[i][k])
                    * (aij[i][i] * aij[k][k]).sqrt();
            }
            t3 /= a_mix;
            t3 -= b[k] / b_mix;

            let t4 = (a_mix / (b_mix * RGASCON * self.t_kel.powf(1.5)))
                * (v_mol / (v_mol + b_mix)).ln();

            let log_phi = t1 + t2 + t3 * t4; // Eq A-8
            self.fug_pi[k] = self.p_bar * log_phi.exp();
        }
    }

    // ---- K-value temperature polynomial (Eq 5) -----------------------
    #[inline]
    fn ktp(&self, k0: f64, v_bar: f64) -> f64 {
        k0 * (self.p_rt0 * v_bar).exp()
    }

    // ---- K_CO2 (Eq 5 + 6) -------------------------------------------
    fn k_co2(&mut self) {
        let coeffs = if self.low_temp {
            if self.co2_sat {
                &[1.169, 1.368e-2, -5.380e-5][..] // Liquid CO2
            } else {
                &[1.189, 1.304e-2, -5.446e-5][..]
            }
        } else {
            &[1.668, 3.992e-3, -1.156e-5, 1.593e-9][..]
        };

        let mut k0 = 10.0_f64.powf(ft(self.deg_c, coeffs));

        if self.scaled {
            let k0_lt = 10.0_f64.powf(ft(self.deg_c, &[1.189, 1.304e-2, -5.446e-5]));
            k0 = self.blended_val(k0_lt, k0);
        }

        self.k_val[0] = self.ktp(k0, self.v_bar[0]);
    }

    // ---- K_H2O (Eq 5 + 6) -------------------------------------------
    fn k_h2o(&mut self) {
        let coeffs = if self.low_temp {
            &[-2.209, 3.097e-2, -1.098e-4, 2.048e-7][..]
        } else {
            &[-2.1077, 2.8127e-2, -8.4298e-5, 1.4969e-7, -1.1812e-10][..]
        };

        let mut k0 = 10.0_f64.powf(ft(self.deg_c, coeffs));

        if self.scaled {
            let k0_lt =
                10.0_f64.powf(ft(self.deg_c, &[-2.209, 3.097e-2, -1.098e-4, 2.048e-7]));
            k0 = self.blended_val(k0_lt, k0);
        }

        self.k_val[1] = self.ktp(k0, self.v_bar[1]);
    }

    // ---- Composition-dependent gammas (Eqs 12, 13, 15) ---------------
    fn calc_gammas(&mut self) {
        if self.low_temp {
            self.gamma = [1.0, 1.0];
            return;
        }
        let am =
            -3.084e-2 * (self.t_kel - 373.15) + 1.927e-5 * (self.t_kel - 373.15).powi(2); // Eq 15
        let x0 = self.x[0];
        self.gamma = [
            (2.0 * am * x0 * (1.0 - x0).powi(2)).exp(),          // Eq 12
            ((am - 2.0 * am * (1.0 - x0)) * x0 * x0).exp(),       // Eq 13
        ];
    }

    // ---- Equilibrium factors A and B' (Eqs 10 and 17) ----------------
    fn a_b(&mut self) {
        let a = self.k_val[1] * self.gamma[1] / self.fug_pi[1]; // Eq 10
        let mut b = self.fug_pi[0]
            / (CONMOLA * self.gamma[0] * self.gamma_prime * self.k_val[0]); // Eq 17
        b = b.clamp(EPS, 1.0 - EPS);
        self.a_eq = a;
        self.b_prime = b;
    }

    // ---- mixMolar helper ---------------------------------------------
    #[inline]
    fn mix_molar(x1: f64, p1: f64, p2: f64) -> f64 {
        x1 * p1 + (1.0 - x1) * p2
    }
}

// =========================================================================
// Halley cubic solver: find ALL real roots of Z^3 + c2*Z^2 + c1*Z + c0 = 0
// =========================================================================
fn halley_all_roots(c2: f64, c1: f64, c0: f64) -> Vec<f64> {
    // Step 1: find one root via Halley iteration starting near inflection
    let mut z = -c2 / 3.0;
    let f0 = z * z * z + c2 * z * z + c1 * z + c0;
    if f0 < 0.0 {
        z += 1.0; // start on vapor side
    }

    let max_iter = 50;
    let tol = 1e-12;
    let mut converged = false;

    for _ in 0..max_iter {
        let f = z * z * z + c2 * z * z + c1 * z + c0;
        let fp = 3.0 * z * z + 2.0 * c2 * z + c1;
        let fpp = 6.0 * z + 2.0 * c2;
        if fp.abs() < 1e-30 {
            break;
        }
        let dz_newton = f / fp;
        let denom = fp - 0.5 * dz_newton * fpp;
        if denom.abs() < 1e-30 {
            break;
        }
        let dz = f / denom;
        z -= dz;
        if dz.abs() < tol {
            converged = true;
            break;
        }
    }

    if !converged {
        // Fallback to Cardano analytical solver
        return cardano_all_roots(c2, c1, c0);
    }

    // Check residual
    let f = z * z * z + c2 * z * z + c1 * z + c0;
    if f.abs() > 1e-6 {
        return cardano_all_roots(c2, c1, c0);
    }

    let r1 = z;

    // Synthetic division to get quadratic Z^2 + e1_q*Z + e0_q
    let e1_q = c2 + r1;
    let e0_q = c1 + r1 * e1_q;

    let mut roots = vec![r1];

    let disc = e1_q * e1_q - 4.0 * e0_q;
    if disc >= 0.0 {
        let sqrt_disc = disc.sqrt();
        let mut r2 = (-e1_q + sqrt_disc) / 2.0;
        let mut r3 = (-e1_q - sqrt_disc) / 2.0;

        // Refine each root with one Halley step
        for r in [&mut r2, &mut r3] {
            let rv = *r;
            let f = rv * rv * rv + c2 * rv * rv + c1 * rv + c0;
            let fp = 3.0 * rv * rv + 2.0 * c2 * rv + c1;
            let fpp = 6.0 * rv + 2.0 * c2;
            if fp.abs() > 1e-30 {
                let dz = f / fp;
                let denom = fp - 0.5 * dz * fpp;
                if denom.abs() > 1e-30 {
                    *r = rv - f / denom;
                }
            }
        }
        roots.push(r2);
        roots.push(r3);
    }

    roots
}

// Cardano analytical solver returning all real roots
fn cardano_all_roots(c2: f64, c1: f64, c0: f64) -> Vec<f64> {
    let p = (3.0 * c1 - c2 * c2) / 3.0;
    let q = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1 + 27.0 * c0) / 27.0;
    let disc = q * q / 4.0 + p * p * p / 27.0;

    if disc < 0.0 {
        // Three real roots (casus irreducibilis)
        let m = 2.0 * (-p / 3.0).sqrt();
        let qpm = (3.0 * q / (p * m)).clamp(-1.0, 1.0);
        let theta = qpm.acos() / 3.0;
        let shift = -c2 / 3.0;

        let r0 = m * theta.cos() + shift;
        let r1 = m * (theta + 4.0 * std::f64::consts::PI / 3.0).cos() + shift;
        let r2 = m * (theta + 2.0 * std::f64::consts::PI / 3.0).cos() + shift;
        vec![r0, r1, r2]
    } else {
        let sqrt_disc = disc.sqrt();
        let pp = -q / 2.0 + sqrt_disc;
        let qq = -q / 2.0 - sqrt_disc;
        let pp_cr = if pp >= 0.0 { pp.cbrt() } else { -(-pp).cbrt() };
        let qq_cr = if qq >= 0.0 { qq.cbrt() } else { -(-qq).cbrt() };
        vec![pp_cr + qq_cr - c2 / 3.0]
    }
}

// =========================================================================
// Main entry point: co2_brine_solubility
//
// Inputs:  p_bar (bar), deg_c (Celsius), ppm (weight NaCl per 1e6 weight brine)
// Returns: (x_co2, y_co2, y_h2o, rho_gas [g/cm3], gas_z)
// =========================================================================
pub fn co2_brine_solubility(
    p_bar: f64,
    deg_c: f64,
    ppm: f64,
) -> (f64, f64, f64, f64, f64) {
    let mut s = SpState::new(p_bar, deg_c, ppm);

    // Clamp pressure floor
    if s.p_bar <= 1.0 {
        s.p_bar = 1.0 + EPS;
    }

    s.t_kel = s.deg_c + CEL2KEL;
    s.p_rt = s.p_bar / (RGASCON * s.t_kel);

    // Reference pressure
    if s.deg_c <= 100.0 {
        s.p0 = 1.0;
    } else {
        s.p0 = ft(
            s.deg_c,
            &[-1.9906e-1, 2.0471e-3, 1.0152e-4, -1.4234e-6, 1.4168e-8],
        );
    }
    s.p_rt0 = (s.p_bar - s.p0) / (RGASCON * s.t_kel);
    s.p_rt = s.p_bar / (RGASCON * s.t_kel);

    // Molality and salt mole fraction
    let fppm = ppm / 1e6;
    s.molal = s.ppm2molality();
    s.x_salt = 2.0 * fppm * MWWAT / (fppm * MWWAT + (1.0 - fppm) * MWSAL);

    // Activity coefficient
    s.gamma_co2();

    // Determine calculation path
    s.calc_type();

    // Component molar volumes
    if s.low_temp {
        s.v_bar = [32.6, 18.1];
    } else {
        let dt = s.t_kel - 373.15;
        s.v_bar = [
            ft(dt, &[32.6, 3.413e-2]),
            ft(dt, &[18.1, 3.137e-2]),
        ];
    }

    // K_H2O first (outside the repeat loop)
    s.k_h2o();

    // The "repeat" loop from the Python code: body always executes exactly once.
    // cubicSolver tracks whether co2_sat state changed (the `repeat` flag), but
    // the Python while-loop structure causes the body to run only once regardless.
    // If co2_sat changed, the flag is noted but does not cause re-execution.
    {
        // K_CO2
        s.k_co2();

        // First estimates of yCO2 and xCO2
        if !s.low_temp {
            let y_co2 = 1.0 - s.est_yh2o();
            let x_co2 = y_co2 / s.k_val[0];
            s.y = [y_co2, 1.0 - y_co2];
            s.x = [x_co2, 1.0 - x_co2];
        } else {
            s.y = [1.0, 0.0];
            s.x = [0.0, 1.0];
        }

        // Mixing rules
        s.a_mix_rk();
        s.b_mix_rk();

        // Solve cubic for molar volume (may update co2_sat state)
        let _repeat = s.molar_volume();
    }

    // Fugacity coefficients * pressure
    s.fug_p();

    // If scaled (99-109 degC): recalculate with low temp coefficients and blend
    if s.scaled {
        let phi_p_ht = s.fug_pi;
        s.low_temp = true;
        s.a_mix_rk();
        s.b_mix_rk();
        // Note: do NOT re-solve cubic -- MolarVolume is not called again
        // The Python code calls self.fugP() which uses self.MolarVol from the
        // high-temp solve above
        s.fug_p();
        s.fug_pi[0] = s.blended_val(s.fug_pi[0], phi_p_ht[0]);
        s.fug_pi[1] = s.blended_val(s.fug_pi[1], phi_p_ht[1]);
        s.low_temp = false;
    }

    s.calc_gammas();
    s.a_b();

    // Calculate results (Eq B-7)
    let denom_yh2o = (1.0 / s.a_eq - s.b_prime) * (2.0 * s.molal + CONMOLA)
        + 2.0 * s.molal * s.b_prime;
    s.y[1] = (1.0 - s.b_prime) * CONMOLA / denom_yh2o;

    s.x[0] = s.b_prime * (1.0 - s.y[1]);
    s.y[0] = 1.0 - s.y[1];

    // Eq B-6
    let m_co2 = s.x[0] * (CONMOLA + 2.0 * s.molal) / (1.0 - s.x[0]);
    // Eq B-3
    s.x_salt = 2.0 * s.molal / (2.0 * s.molal + CONMOLA + m_co2);
    s.x[1] = (1.0 - s.x[0] - s.x_salt).clamp(0.0, 1.0);

    // Iterative refinement for high-temp path
    if !s.low_temp {
        let mut err = 1.0;
        let mut iter_num = 0;
        while err > EPS && iter_num < 100 {
            let yh2o_last = s.y[1].max(EPS);

            // Mixing rules for updated compositions
            s.a_mix_rk();
            s.b_mix_rk();

            // Fugacity
            s.fug_p();

            if s.scaled {
                let phi_p_ht = s.fug_pi;
                s.low_temp = true;
                s.a_mix_rk();
                s.b_mix_rk();
                s.fug_p();
                s.fug_pi[0] = s.blended_val(s.fug_pi[0], phi_p_ht[0]);
                s.fug_pi[1] = s.blended_val(s.fug_pi[1], phi_p_ht[1]);
                s.low_temp = false;
            }

            s.calc_gammas();
            s.a_b();

            let denom = (1.0 / s.a_eq - s.b_prime) * (2.0 * s.molal + CONMOLA)
                + 2.0 * s.molal * s.b_prime;
            s.y[1] = ((1.0 - s.b_prime) * CONMOLA / denom).clamp(EPS, 1.0 - EPS);
            s.x[0] = (s.b_prime * (1.0 - s.y[1])).clamp(EPS, 1.0 - EPS);
            s.y[0] = 1.0 - s.y[1];

            let m_co2 = s.x[0] * (CONMOLA + 2.0 * s.molal) / (1.0 - s.x[0]);
            s.x_salt = 2.0 * s.molal / (2.0 * s.molal + CONMOLA + m_co2);
            s.x[1] = 1.0 - s.x[0] - s.x_salt;

            err = (s.y[1] / yh2o_last - 1.0).abs();
            iter_num += 1;
        }
    }

    // Re-compute gas phase density
    let mw_gas = SpState::mix_molar(s.y[0], MWCO2, MWWAT);
    let rho_gas = mw_gas / s.molar_vol;
    let gas_z = s.molar_vol * s.p_rt;

    (s.x[0], s.y[0], s.y[1], rho_gas, gas_z)
}

// =========================================================================
// Unit tests
// =========================================================================
#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol || (a - b).abs() / b.abs().max(1e-30) < tol
    }

    #[test]
    fn test_ft_polynomial() {
        // FT([1.0, 2.0, 3.0], t=2) = 1 + 2*2 + 3*4 = 17
        assert!((ft(2.0, &[1.0, 2.0, 3.0]) - 17.0).abs() < 1e-12);
    }

    #[test]
    fn test_halley_single_root() {
        // Z^3 - 6Z^2 + 11Z - 6 = 0  => roots 1, 2, 3
        let roots = halley_all_roots(-6.0, 11.0, -6.0);
        assert!(roots.len() == 3);
        let mut sorted = roots.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(approx_eq(sorted[0], 1.0, 1e-8));
        assert!(approx_eq(sorted[1], 2.0, 1e-8));
        assert!(approx_eq(sorted[2], 3.0, 1e-8));
    }

    #[test]
    fn test_low_temp_basic() {
        // 100 bar, 50 degC, 0 ppm -- should produce reasonable xCO2
        let (x_co2, y_co2, y_h2o, rho_gas, gas_z) = co2_brine_solubility(100.0, 50.0, 0.0);
        assert!(x_co2 > 0.0 && x_co2 < 0.1, "xCO2 out of range: {}", x_co2);
        assert!(y_co2 > 0.9 && y_co2 <= 1.0, "yCO2 out of range: {}", y_co2);
        assert!(y_h2o >= 0.0 && y_h2o < 0.1, "yH2O out of range: {}", y_h2o);
        assert!(rho_gas > 0.0, "rhoGas must be positive: {}", rho_gas);
        assert!(gas_z > 0.0 && gas_z < 2.0, "gasZ out of range: {}", gas_z);
    }

    #[test]
    fn test_high_temp() {
        // 200 bar, 150 degC, 30000 ppm
        let (x_co2, y_co2, y_h2o, rho_gas, gas_z) = co2_brine_solubility(200.0, 150.0, 30000.0);
        assert!(x_co2 > 0.0 && x_co2 < 0.1, "xCO2 out of range: {}", x_co2);
        assert!(y_co2 > 0.5 && y_co2 <= 1.0, "yCO2 out of range: {}", y_co2);
        assert!(y_h2o > 0.0, "yH2O should be > 0 at high temp: {}", y_h2o);
        assert!(rho_gas > 0.0, "rhoGas must be positive: {}", rho_gas);
        assert!(gas_z > 0.0, "gasZ must be positive: {}", gas_z);
    }

    #[test]
    fn test_scaled_range() {
        // 150 bar, 105 degC, 10000 ppm -- falls in blended range
        let (x_co2, y_co2, y_h2o, rho_gas, gas_z) = co2_brine_solubility(150.0, 105.0, 10000.0);
        assert!(x_co2 > 0.0 && x_co2 < 0.1, "xCO2 out of range: {}", x_co2);
        assert!(y_co2 > 0.5, "yCO2 should be dominant: {}", y_co2);
        assert!(y_h2o >= 0.0, "yH2O must be non-negative: {}", y_h2o);
        assert!(rho_gas > 0.0, "rhoGas must be positive: {}", rho_gas);
        assert!(gas_z > 0.0, "gasZ must be positive: {}", gas_z);
    }
}
