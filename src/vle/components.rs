/// Component properties database for Soreide-Whitson VLE engine.
///
/// 16 components: H2O + 15 gas species.
/// Critical properties from S&W 1992 Table 5; H2 from this work.

/// Physical constants
pub const R_GAS: f64 = 8.314462; // J/(mol·K)
pub const OMEGA_A: f64 = 0.45724;
pub const OMEGA_B: f64 = 0.07780;
pub const MW_NACL: f64 = 58.44;
pub const MW_H2O: f64 = 18.015;

/// H2 critical temperature for BIP correlations (NIST)
pub const BIP_TC_H2: f64 = 33.145; // K

/// Component critical properties.
#[derive(Debug, Clone, Copy)]
pub struct ComponentProperties {
    pub tc: f64,    // Critical temperature (K)
    pub pc: f64,    // Critical pressure (Pa)
    pub omega: f64, // Acentric factor
    pub tb: f64,    // Normal boiling point (K) - for Sechenov
    pub mw: f64,    // Molecular weight (g/mol)
}

/// Component name constants (indices into COMPONENT_NAMES / COMPONENT_DB).
pub const IDX_H2O: usize = 0;
pub const IDX_H2: usize = 1;
pub const IDX_CO2: usize = 2;
pub const IDX_H2S: usize = 3;
pub const IDX_N2: usize = 4;
pub const IDX_CH4: usize = 5;
pub const IDX_C2H6: usize = 6;
pub const IDX_C3H8: usize = 7;
pub const IDX_IC4H10: usize = 8;
pub const IDX_NC4H10: usize = 9;
pub const IDX_IC5H12: usize = 10;
pub const IDX_NC5H12: usize = 11;
pub const IDX_NC6H14: usize = 12;
pub const IDX_NC7H16: usize = 13;
pub const IDX_NC8H18: usize = 14;
pub const IDX_NC10H22: usize = 15;

pub const NUM_COMPONENTS: usize = 16;

/// Ordered component names matching COMPONENT_DB indices.
pub const COMPONENT_NAMES: [&str; NUM_COMPONENTS] = [
    "H2O", "H2", "CO2", "H2S", "N2", "CH4",
    "C2H6", "C3H8", "iC4H10", "nC4H10",
    "iC5H12", "nC5H12", "nC6H14", "nC7H16", "nC8H18", "nC10H22",
];

/// Component database: Tc, Pc, omega, Tb, MW
/// Order matches COMPONENT_NAMES.
pub const COMPONENT_DB: [ComponentProperties; NUM_COMPONENTS] = [
    // H2O
    ComponentProperties { tc: 647.3, pc: 22.12e6, omega: 0.3434, tb: 373.15, mw: 18.015 },
    // H2
    ComponentProperties { tc: 33.145, pc: 1.2964e6, omega: -0.219, tb: 20.3, mw: 2.016 },
    // CO2
    ComponentProperties { tc: 304.2, pc: 7.38e6, omega: 0.2273, tb: 194.7, mw: 44.01 },
    // H2S
    ComponentProperties { tc: 373.2, pc: 8.94e6, omega: 0.1081, tb: 212.8, mw: 34.082 },
    // N2
    ComponentProperties { tc: 126.1, pc: 3.40e6, omega: 0.0403, tb: 77.36, mw: 28.014 },
    // CH4
    ComponentProperties { tc: 190.6, pc: 4.60e6, omega: 0.0108, tb: 111.66, mw: 16.043 },
    // C2H6
    ComponentProperties { tc: 305.4, pc: 4.88e6, omega: 0.0986, tb: 184.6, mw: 30.07 },
    // C3H8
    ComponentProperties { tc: 369.8, pc: 4.25e6, omega: 0.1524, tb: 231.1, mw: 44.097 },
    // iC4H10
    ComponentProperties { tc: 408.1, pc: 3.65e6, omega: 0.1770, tb: 261.4, mw: 58.123 },
    // nC4H10
    ComponentProperties { tc: 425.2, pc: 3.80e6, omega: 0.1931, tb: 272.7, mw: 58.123 },
    // iC5H12
    ComponentProperties { tc: 460.4, pc: 3.38e6, omega: 0.2270, tb: 301.0, mw: 72.15 },
    // nC5H12
    ComponentProperties { tc: 469.6, pc: 3.37e6, omega: 0.2510, tb: 309.2, mw: 72.15 },
    // nC6H14
    ComponentProperties { tc: 507.4, pc: 3.01e6, omega: 0.2990, tb: 341.9, mw: 86.18 },
    // nC7H16
    ComponentProperties { tc: 540.3, pc: 2.74e6, omega: 0.3490, tb: 371.6, mw: 100.2 },
    // nC8H18
    ComponentProperties { tc: 568.8, pc: 2.49e6, omega: 0.3980, tb: 398.8, mw: 114.2 },
    // nC10H22
    ComponentProperties { tc: 617.7, pc: 2.10e6, omega: 0.4900, tb: 447.3, mw: 142.3 },
];

/// Look up component index by name. Returns None if not found.
pub fn comp_index(name: &str) -> Option<usize> {
    COMPONENT_NAMES.iter().position(|&n| n == name)
}

/// Build arrays of Tc, Pc, omega, Tb for a given set of component indices.
pub struct ComponentArrays {
    pub tc: Vec<f64>,
    pub pc: Vec<f64>,
    pub omega: Vec<f64>,
    pub tb: Vec<f64>,
    pub mw: Vec<f64>,
}

impl ComponentArrays {
    pub fn from_indices(indices: &[usize]) -> Self {
        let n = indices.len();
        let mut tc = Vec::with_capacity(n);
        let mut pc = Vec::with_capacity(n);
        let mut omega = Vec::with_capacity(n);
        let mut tb = Vec::with_capacity(n);
        let mut mw = Vec::with_capacity(n);
        for &idx in indices {
            let c = &COMPONENT_DB[idx];
            tc.push(c.tc);
            pc.push(c.pc);
            omega.push(c.omega);
            tb.push(c.tb);
            mw.push(c.mw);
        }
        ComponentArrays { tc, pc, omega, tb, mw }
    }
}
