# Correctness — brine + plyasunov

## Critical (bug, user gets wrong answer)

- **brine/_lib_vle_engine.py:735-753 vs src/vle/flash.rs:209-221** — Python/Rust divergence in Sechenov `ks` for CO2 and H2S. Python `get_sechenov_ks()` routes CO2/H2S to Dubessy/Akinfiev specialized models with empirical shifts; Rust only provides S&W Eq 8 fallback. With Rust accelerator + `framework='proposed'`, ks diverges 0.03–0.05 abs, K-values shift 5–10%, underestimating dissolved gas solubility. Silent.
  **Fix:** disable Rust for `proposed`, or port Dubessy/Akinfiev, or warn on fallback.

- **brine/brine.py:1267-1273 (garciaDensity)** — Garcia Eq. 18 computes `xRat = xCO2/(1-xCO2) → ∞` as `xCO2 → 1`; no guard. Prior 2026-04-12 review noted this guard should exist; it is **not present in current code**. Silent NaN/zero density at high CO2 fractions (CO2 storage scenarios).
  **Fix:** `if xCO2 > 1 - 1e-6: warn + return unsaturated brine density`. Also guard salt_ratio denom near-zero at line 1241.

- **plyasunov/plyasunov_model.py:254-273** — `_compute_ai` uses `theta**last_exp` with `last_exp ∈ {-9,...,-6}`. For T < 50 K, `theta < 0.077`, giving `theta^-9 > 5e14`. No bounds checking — silently extrapolates beyond IAPWS-IF97 Region 1 validity (273–623 K).
  **Fix:** document valid range, clamp-with-warning at boundaries.

## High (edge case, silent failure)

- **brine/_lib_vle_engine.py:127-175 (alpha functions)** — `alpha_water_soreide` clamps Tr to min 0.01 but does NOT prevent Tr > 1. For T > Tc_water (647.3 K), sqrt_alpha^2 blows up. `alpha_water_mc3` cubic can go negative. Used in Spycher-Pruess CO2 mixture at T > 100 C — fugacity diverges above ~130 C.
  **Fix:** explicit Tr range check: `if Tr < 0.01 or Tr > 1.2: warn + clamp`.

- **brine/_lib_vle_engine.py:2143-2176 (flash_tp), 1110-1142 (solve_rachford_rice)** — RR solver clips K-values AFTER the single-phase decision. Extreme initial K values can fool single-phase detection (two-phase false positive), and compositions come back non-normalized (1e-15 noise). Also line 1041 perturbs K=1 to 1+1e-12, causing component to vanish from RR eq.
  **Fix:** normalize K before single-phase check.

- **brine/brine.py:1127-1133 (CO2_Brine_Mixture Spycher iter)** — Non-convergence warns but does NOT return/store a `converged` flag. Downstream `brine_props()` uses the non-converged (x,y) silently at lines 1146-1149. `SoreideWhitson.flash_tp` returns a `converged` bool — `CO2_Brine_Mixture` does not.
  **Fix:** expose `.converged` attribute or return flag.

## Medium (validity / robustness)

- **SoreideWhitson framework × salinity_method combinations** — orthogonal params with no validation: `framework='proposed' + salinity_method='embedded'` silently falls back to gamma_phi; `framework='sw_original' + salinity_method='gamma_phi'` double-counts salinity (embedded in kij AND in gamma). **Fix:** validate combinations in `calc_equilibrium()`.

- **IAPWS-IF97 above 100 MPa** — `_rho_if97()` called from `brine.py:185, 223, 1230` with no upper-pressure check. Valid range is Psat(T) ≤ P ≤ 100 MPa. Silent extrapolation for deep CO2/H2 storage. **Fix:** document + optional clamp/warn.

- **Sechenov ks sign convention** — docstring says `ks>0 ⟹ salting-out`; Dubessy CO2 `ks` can flip sign across T/m range, Akinfiev H2S extracted by difference. No assertion that `gamma >= 0.1`. Sign flip causes 10x K-value error.

- **_lib_salting_library.py:215-266 (Akinfiev H2S tables)** — `np.interp` default flat extrapolation on T=298.15–523.15 K table; T outside range silently clamps. Table at T=298.15 has NaN (missing data).

## Low (nitpick)

- **brine/brine.py:59-100 (*_ARR arrays)** — 1-based indexing with `[0]` placeholder (McCain Eq 4.1 compatibility) is undocumented. Add comment.
- **plyasunov P_IN_MAP** — gas name normalization uppercase; underscore/hyphen variants silently fail. Add docstring.
- **tests/test_brine.py** — no tests for xCO2>0.99 Garcia, T>Tc_water alpha, invalid framework/salinity combos, P>100 MPa, Spycher non-convergence behavior.

## Top 3 priorities

1. **Garcia density singularity** — silent NaN at high xCO2 is a silent-failure regression (prior review noted fix was in place).
2. **Python/Rust ks divergence** — any user running with Rust accelerator on `proposed` framework gets wrong answers silently.
3. **Spycher non-convergence not exposed** — users have no way to detect it programmatically.
