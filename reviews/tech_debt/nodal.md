# Tech Debt — nodal

Scope: `pyrestoolbox/nodal/nodal.py` (2,132 lines, single file) + `src/vlp/*.rs` (9 files, ~4,100 lines total). `src/nodal/` does not exist — all Rust acceleration lives under `src/vlp/`.

## High

### H1. Rust side has no segment-march scaffold — 8 near-duplicate segment loops
`src/vlp/segment_gas.rs` defines `hb_fbhp_gas`, `wg_fbhp_gas`, `gray_fbhp_gas`, `bb_fbhp_gas` (4 functions, 545 lines). `src/vlp/segment_oil.rs` mirrors this with 4 oil variants (593 lines). Each has its own full preamble (Sutton Tc/Pc, mass check, low-rate fallback, geometry, `calc_segments`, 2-iter pressure loop, per-step PVT, IFT, velocities). The Python `_segment_march_gas` / `_segment_march_oil` scaffold does NOT have a Rust analogue. ~800 lines of copy-paste across the 8 functions — the most recent change (HB temp fix) had to be propagated to both Python and 2 Rust files. The decorator `@rust_accelerated` on all 8 Python wrappers masks the Rust duplication at the call site but doesn't fix it. Refactor `src/vlp/segment_gas.rs` / `segment_oil.rs` to match the Python scaffold (single march function + 4 closure-based gradient callbacks, or enum-dispatched gradient fn pointer).

### H2. Python-only constants pass — Rust is still magic-number heavy
The Tier 3 "named constants with paper citations" pass was Python-only. Rust `src/vlp/` retains `0.0765`, `62.4`, `5.615`, `86400`, `141.5`, `131.5`, `0.4846`, `0.5351`, `0.5824`, `0.0868`, `0.0173`, `0.0609`, `3.7680`, `3.5390`, `1.6140`, `2.960`, `0.3050`, `-0.4473`, `0.0978`, `0.924` (Payne), `316.0`, `0.0009252`, `-2.4684`, `0.10`, `-1.4516`, `0.5`, `-6.738`, `0.9116257`, `-4.821756`, `1232.25`, `-22253.58`, `116174.3`, `1.938`, `120.872`, `0.15726`, `1.22`, `2.9`, `28.97`, `1.11591`, `-0.00305`, `38.085`, `-0.259`, `32.0436`, `-1.1367`, `75.0`, `-1.108`, `0.349`, etc., all inlined. `src/vlp/holdup_bb.rs`, `holdup_hb.rs`, `ift.rs`, `segment_gas.rs`, `segment_oil.rs` need a dedicated `constants.rs` or per-file `const` block mirroring Python's named constants.

### H3. Rust parity test suite is threadbare (8 tests, 1 flow scenario, loose rtol)
`test_rust_acceleration.py::TestNodalVLPEquivalence` covers exactly one gas scenario (`thp=500, qg_mmscfd=5, sg=0.75, 10000-ft vertical, T=100-200`) and one oil scenario, across 4 methods × 2 phases = 8 tests, all at `RTOL_LOOSE` with `atol=1.0`. This gives no coverage of:
- High-deviation / horizontal wells (theta != pi/2)
- Multi-segment completions
- Static-column fallback (low rate)
- Condensate dropout (pr>0)
- Injection path (`injection=True`)
- Metric unit path
- BB TRANSITION flow pattern (correctness review flagged a `** 4` bug that exists in both Python AND Rust — this parity test set would not catch divergence because both are wrong together)

Parity tests should be expanded to a parametric sweep. Related: `_compare_fbhp_gas` uses `atol=1.0 psi` — large BHP differences fit inside this band.

## Medium

### M1. Residual duplication across 4 gradient callbacks
Each of `_hb_gradient_gas`, `_wg_gradient_gas`, `_gray_gradient_gas`, `_bb_gradient_gas` recomputes the same blocks independently:
- `rho_s = rho_l * hl + rho_g * (1-hl)` (HB, GRAY, BB — 3 copies)
- `mu_ns = mu_l * lambda_l + mu_g * (1-lambda_l)` then `mu_ns_lbfts = mu_ns * _CP_TO_LBFTS` (GRAY, BB)
- `n_re = rho_ns * v_m * diam_ft / mu_ns_lbfts` (GRAY, BB — 2 copies, identical)
- `eps_d = rough_ft / diam_ft` (BB line 1589; HB uses `rough / tid` on line 978 — inconsistent units for identical intent)
- `ek = gm*v_sg / (GC*p*144)` (GRAY) vs `ek = rho_s*v_m*v_sg / (GC*p*144)` (BB) — two similar acceleration-term denominators, boilerplate clamped `denom = max(1-ek, 0.1)`
Extract to helpers in the shared-march layer (`_mix_density`, `_no_slip_viscosity_lbfts`, `_reynolds`, `_accel_denom`).

### M2. BB TRANSITION branch recomputes `_bb_horizontal_holdup` redundantly (nodal.py:1555–1575)
Lines 1555–1556 compute `hl0_seg` and `hl0_int`. Then lines 1571 and 1574 pass `_bb_horizontal_holdup(..., _SEGREGATED) * _BB_PAYNE` and `_bb_horizontal_holdup(..., _INTERMITTENT) * _BB_PAYNE` again into `_bb_inclination_correction` — two redundant re-evaluations of the same polynomial. The same duplication exists in `src/vlp/segment_gas.rs:495, 499`. Pass the already-computed `hl0_seg` / `hl0_int` through. (Flagged as "duplicative but correct" in correctness/nodal.md M3; re-listing here as tech debt since the fix is trivial.)

### M3. State-dict payload has ~11 keys never consumed by any callback
`_segment_march_gas` and `_segment_march_oil` populate `s['qg_mmscfd']`, `s['osg']`, `s['gsg']`, `s['qw_bwpd']`, `s['oil_vis_loc']`, `s['temp_r']`, `s['zee']`, `s['qo_loc']`, `s['rs_est']`, `s['tc']`, `s['pc']`, `s['api']`. Grepping the 4 callbacks shows none of these are read. Dict allocations run once per segment × 2 pressure iters × ~100 segments × N rates → modest GC pressure. Strip unused keys, or switch to a `dataclass` / named tuple for the state container (saves dict hashing overhead too). If Rust ever grows a callback-based scaffold, an explicit struct is mandatory — good to align now.

### M4. Duplicated metric-boundary and oil_pvt-extraction blocks across 3 public APIs
- Metric conversion preamble at `fbhp:1687-1705`, `outflow_curve:1795-1810`, `operating_point:1995-2008` — 3 near-identical blocks (thp→psia, pr→psia, cgr/qw/qt/gor/pb/rsb conversions).
- `oil_pvt` attribute extraction at `fbhp:1713-1721`, `outflow_curve:1812-1818`, `operating_point:2013-2019` — 3 copies of the same 7-line block.
- Per-well-type rate conversion scaffolding inside `outflow_curve` (lines 1822–1828, 1837, 1844, 1850) duplicates the same gas-vs-oil unit logic as `_convert_vlp_to_metric` (which itself exists as a separate helper at 2115).

Extract to `_convert_vlp_inputs_to_oilfield(...)`, `_extract_oil_pvt_params(oil_pvt)` helpers.

### M5. Copied PVT correlations with their own magic numbers
`_z_factor`, `_gas_viscosity`, `_gas_density`, `_water_viscosity`, `_standing_rs`, `_velarde_rs`, `_oil_density_mccain`, `_oil_viscosity_full` (lines 608-760) are performance-motivated copies of logic that also exists in `pyrestoolbox/gas/` and `pyrestoolbox/oil/`. The nodal copies retain all their magic numbers inline (Hall-Yarborough `-0.06125, 14.76, 9.76, 4.58, 90.7, 242.2, 42.4, 2.18, 2.82`; LGE `9.4, 0.02, 209.0, 19.0, 3.5, 986.0, 2.4, 0.2`; McCain `-49.893, 85.0149, 3.70373, …, 4600.0, 73.71`; Standing `0.00091, 0.0125, 18.2, 1.2048`; Velarde 18 coefficients; Beggs-Robinson `3.0324, 0.02023, 1.163, 10.715, 0.515, 5.44, 0.338`). The main gas/oil modules extracted these; the nodal copies did not. Either reference the canonical module (if the perf hit is tolerable — benchmark first) or duplicate the named constants too.

### M6. Condensate dropout model has no formalized deprecation / warning
`_condensate_dropout` (nodal.py:842) is the known-naive linear model flagged as Tier 4 R&D (prototype moved to `/home/mark/projects/CondensateDropout` per MEMORY.md). No docstring warning, no user-facing note in `fbhp` docstring beyond `"0 disables"`, no `DeprecationWarning`, no TODO in source. Either add a documentation note acknowledging it's a simple linear screening model, or gate it behind a kwarg (`condensate_model='linear'` default) so future replacement is non-breaking.

### M7. Unused params in `fbhp` / `outflow_curve` / `operating_point` are silently accepted
- `gas_pvt` parameter present on `fbhp`, `outflow_curve`, `operating_point` — docstring admits "unused by VLP methods directly, reserved for future use" (nodal.py:1684). Reserved unused params are dead API. Either wire it through (composition passed to the internal `_z_factor`) or drop until used.
- `_hb_fbhp_oil`, `_wg_fbhp_oil`, `_gray_fbhp_oil`, `_bb_core_oil` accept both `rsb_scale` and `rsb_frac` as separate kwargs that ultimately collapse to the same physical quantity inside `_segment_march_oil` (line 1295 passes `rsb_scale`, but line 1736 passes `rsb_scale=rsb_frac`). Naming is inconsistent and the coupling is confusing.

## Low

### L1. Zero type hints across the entire module
2,132 lines, 52 functions, 4 public APIs, 3 data classes — **not a single `-> type` return annotation or parameter hint**. (`Grep "-> "` returns only comment-style unit annotations, not type hints.) For an LLM-agent-facing library (per CLAUDE.md Overview), type hints materially help downstream IDE inference and static analysis. At minimum: annotate the public API (`fbhp -> float`, `outflow_curve -> NodalResult`, `ipr_curve -> NodalResult`, `operating_point -> NodalResult`, and the 3 data class `__init__` signatures).

### L2. Unused imports
`from pyrestoolbox.classes import vlp_method, class_dic` (line 49) — neither is referenced anywhere in the file (grep confirms). `validate_methods` returns the already-resolved Enum.
`from pyrestoolbox.constants import` (lines 52-60) imports but never uses: `degf_to_degc`, `SCF_PER_STB_TO_SM3_PER_SM3`, `SM3_PER_SM3_TO_STB_PER_MSCF`, `STB_PER_MSCF_TO_SM3_PER_SM3`, `STB_PER_MMSCF_TO_SM3_PER_SM3`, `M3_TO_BBL`, `BBL_TO_M3`, `SM3_TO_MSCF`, `D_PER_MSCF_TO_D_PER_SM3`. 9 unused constants plus 2 unused classes imports.

### L3. Unused internal helper
`_gas_density` (line 671) is defined but never called — the same computation is inlined as `_MW_AIR * gsg * p_avg / (zee * _R_GAS * temp_r)` in both `_segment_march_gas:1245` and `_segment_march_oil:1345`. Either call the helper from both segment marches (DRY) or delete it.

### L4. Duplicate `_LBFT3_TO_KGM3` constant definition
Defined twice (line 72 and line 104) with the same value (`16.01846`). Same for `_IN_TO_M` (line 103, matches imported `IN_TO_MM`/12 geometry but is a separate literal). Consolidate the constant block near the top.

### L5. File size — single-file module at 2,132 lines
Oil module hit this cap and was split into 10 sub-files under `pyrestoolbox/oil/`. Nodal is the largest remaining single-file module. Natural split: `_pvt_helpers.py` (lines 608-760), `_ift.py` (767-806), `_friction.py` (813-825), `_gradients.py` (all 4 callbacks), `_segment_march.py`, `nodal.py` (data classes + public API). Not urgent — the file is logically segmented with banner comments — but consistency with the oil sub-file pattern is worth noting.

### L6. Magic number `_GAS_COL_K = 0.01875`
Line 122, used at line 1355 in the oil segment march's fallback degenerate-mass branch. No paper citation, no derivation comment. Likely `g*MW_air / (R * 144)` in some unit combination, but should be documented or derived inline.

### L7. `oil_pvt.rsb_frac`, `oil_pvt.vis_frac` extracted without fallback
Line 1718-1719: `vis_frac = oil_pvt.vis_frac`, `rsb_frac = oil_pvt.rsb_frac`. If a future OilPVT drops those attributes, this silently errors. Use `getattr(oil_pvt, 'vis_frac', 1.0)`. Minor defensive-coding nit.

### L8. `operating_point` silent bisection-failure returns `op_rate=0` (also correctness M)
Tech-debt angle: silently setting `op_rate=0.0` on `bisect_solve` failure (line 2073) wrecks diagnosability. At minimum log the exception via `warnings.warn`. (Correctness review already has this under High, so just a cross-ref.)

---

## Summary

**22 items: 3 High, 7 Medium, 8 Low.** Core finding: the Python-side refactor (Tier 3 constants, `_segment_march` scaffold, `@rust_accelerated` decorator) was not propagated to Rust, leaving Rust with 8 copy-pasted segment loops and inline magic numbers. Combined with a thin parity test suite (8 tests, single scenario, loose tolerance), future algorithm changes risk silent Python/Rust divergence.

**Top 3 priorities:**
1. **H1** — Factor Rust `src/vlp/segment_gas.rs` + `segment_oil.rs` into a shared scaffold mirroring Python; ~800 lines of duplication is the single biggest maintenance liability.
2. **H3** — Expand `TestNodalVLPEquivalence` to a parametric sweep (deviation, multi-segment, static-column fallback, metric, injection, condensate dropout). Current suite would not catch the BB `** 4` bug because Python and Rust are wrong together.
3. **H2 / M5** — Rust constants extraction to match Python's named-constants pass; deduplicate PVT-helper magic numbers (Python and Rust both).
