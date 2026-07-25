# Tech Debt — gas

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/gas/gas.py` (2195 lines, 100 KB, one monolithic file). Prior tech-debt pass (2026-04-12) introduced the named-constants block, `@rust_accelerated` decorator, `_metric_to_field_pvt` helper, and `_prepare_gas_rate_inputs`/`_compute_delta_mp` refactors; items already cleared there are not re-flagged. Nothing overlaps with `reviews/correctness/gas.md` or `reviews/api_ux/gas.md`.

## High (significant maintainability risk / frequent-touch area)

### H-1 — `gas.py` is a 2195-line monolith spanning 4 distinct concerns
Oil was already split into 10 sub-files under `oil/` per the prior Tier 3 pass — gas was not. Four clearly separable concerns in one file:
1. **PVT core** (gas_z, gas_ug, gas_bg, gas_cg, gas_den, gas_tc_pc, gas_sg, gas_dmp, gas_ponz2p, gas_grad2sg, gas_fws_sg, gas_water_content, GasPVT) — lines ~80–1712
2. **BNS/PR-EOS internals** (_BNS_*, _BIP_*, _calc_bips_fast, _cardano_cubic, _halley_cubic_vec) — lines ~573–702 (could live in `_bns_eos.py`)
3. **Darcy / rate engine** (darcy_gas, gas_rate_radial, gas_rate_linear, _compute_delta_mp, _prepare_gas_rate_inputs) — lines ~193–433
4. **Hydrate prediction** (HydrateResult, gas_hydrate, _motiee_hft, _towler_mokhatab_hft, _hydrate_formation_press, _ostergaard_depression, _required_concentration, and 5 dicts) — lines ~1715–2195 (nearly 500 lines; fully self-contained and the obvious first candidate to split into `_hydrate.py`)

The hydrate section is the lowest-risk carve-out: zero cross-calls into PVT except `gas_water_content`/`SoreideWhitson`, its own constants (`_PSI_TO_KPA`, `_HAMMERSCHMIDT`, `_OSTERGAARD`, `_INHIBITOR_DENSITY`, `_GCM3_TO_LB_PER_GAL`, `_WATER_LB_PER_STB`, `_MAX_WT_PCT`, `_MOTIEE`, `_TOWLER`, `_HFP_P_LO/HI`), and only ~4 public entry points (`gas_hydrate`, `HydrateResult` + re-exported enums).

### H-2 — Copy-paste method-resolution block duplicated 10× across public functions
The sequence
```python
zmethod, cmethod = _h2_method_override(h2, zmethod, cmethod)
zmethod, cmethod = validate_methods(["zmethod", "cmethod"], [zmethod, cmethod])
```
appears at lines 748/754, 1036/1042, 1165/1166, 1227/1228, 1277/1278, 1349/1350, 1465/1466, 1526/1527, and the GasPVT ctor (1634/1635). This is the exact pattern `_prepare_gas_rate_inputs` (line 229) was introduced to solve — but the fix was applied only to the two rate functions. A single `_resolve_methods(h2, zmethod, cmethod)` helper (returning the validated Enum pair) would eliminate ~20 redundant lines and make adding a new Z-method or composition override a one-edit change instead of ten.

Worse, `gas_z` (749-752) and `gas_ug` (1037-1040) *also* contain a BNS cross-coupling block:
```python
if cmethod == 'BNS':
    zmethod = 'BNS'
elif zmethod == 'BNS':
    cmethod = 'BNS'
```
…which is *missing* from `gas_bg`, `gas_cg`, `gas_den`, `gas_ponz2p`, `gas_dmp`, `gas_grad2sg`. That asymmetry is either a latent bug or unnecessary in gas_z/gas_ug — but at minimum it is an inconsistency that a single shared helper would have prevented.

### H-3 — `_metric_to_field_pvt` is applied inconsistently; metric conversion is duplicated 6×
Noted as Low in api_ux/gas.md L-6, but from a tech-debt standpoint this is a more substantive copy-paste problem. `gas_z`, `gas_bg`, `gas_cg`, `gas_den` use the helper (one-liner). `gas_ug` (1025-1031), `gas_ponz2p` (1341-1347), `gas_grad2sg` (1442-1449), `gas_dmp` (1511-1518), `gas_rate_radial` (297-310), `gas_rate_linear` (381-391) each inline the same 5-7 line conversion block. Two reasons for the divergence:
- `gas_ug` has extra parameters (`zee`) that aren't in the helper — easy to extend the helper
- rate functions need `r_w`/`length`/`area`/`D` conversions the helper doesn't cover — add a second `_metric_to_field_rate` helper or a kwargs form

The recurring idiom `p = np.asarray(p) * BAR_TO_PSI if not isinstance(p, (int, float)) else p * BAR_TO_PSI` (lines 182, 298, 299, 382, 383, 1026, 1342, 1643) is also over-defensive — `np.asarray(p) * BAR_TO_PSI` works for scalars and arrays alike.

## Medium (worth scheduling a cleanup)

### M-1 — Unused imports at module level
- Line 64 `from typing import Tuple` — used only at line 445 for a return annotation; Python 3.9+ allows lowercase `tuple` (already noted as L-5 in api_ux; leaving for cp38 compat is defensible, listing here for completeness of dead-code scan)
- Line 67 `import pandas as pd` — `pd` has **zero** references in the file. Dead import.
- Line 68 — unused classes imported: `pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic` (9 names). Only `z_method, c_method, hyd_method, inhibitor` are referenced.
- Line 71-79 unused constants: `tscr, FT_TO_M, MM_TO_IN, IN_TO_MM, SQFT_TO_SQM, KGM3_TO_LBCUFT, INVBAR_TO_INVPSI, BAR2CP_TO_PSI2CP, PSIFT_TO_BARM, D_PER_MSCF_TO_D_PER_SM3, SM3_PER_SM3_TO_STB_PER_MSCF, SM3_TO_MSCF` (12 names). Each appears exactly once — on the import line itself.

Net: 22 unused names on the import lines. Cleanable in a single commit.

### M-2 — Duplicated MW magic numbers despite named constants already imported
- Line 500 (SUT branch in `gas_tc_pc`): `sg_hc = (sg - (n2 * 28.01 + co2 * 44.01 + h2s * 34.1) / MW_AIR) / hc_frac`. Uses literals `28.01, 44.01, 34.1` when `MW_N2, MW_CO2, MW_H2S` (= `28.014, 44.01, 34.082`) are imported from constants. Also the H2S value `34.1` does not match `MW_H2S = 34.082` — a genuine minor inconsistency. (Not a correctness issue worth re-flagging; SUT's original paper used these rounded values. Comment should cite the paper value so a later reader doesn't "helpfully" harmonise.)
- Lines 1286-1292 (`gas_den` pure-component branch): `if co2 == 1: m = 44.01 / if h2s == 1: m = 34.082 / if n2 == 1: m = 28.014 / if h2 == 1: m = 2.016`. Same four MW literals duplicated; already in `MW_CO2, MW_H2S, MW_N2, MW_H2` constants. Also uses chained `if` (not `elif`) — cheap correctness-adjacent inefficiency but not a bug since mole fractions sum to 1.
- `_BNS_MWS = np.array([44.01, 34.082, 28.014, 2.016, 0])` (line 574) also duplicates the same four values. Could be built from the imported constants: `np.array([MW_CO2, MW_H2S, MW_N2, MW_H2, 0.0])`.

### M-3 — `tc * 1.8` K→degR conversion unnamed and repeated 8×
Lines 185, 308, 389, 464, 1029, 1345, 1447, 1516 all use bare `1.8` with a `# K to deg R` comment. Either add `K_TO_DEGR = 9.0/5.0` to `pyrestoolbox.constants` (matches the `degc_to_degf` pattern) or a local `_K_TO_DEGR` constant. Pure nit-level magic number, but 8 occurrences.

### M-4 — `sg_hc` hydrocarbon-SG calculation is computed 3 different ways in 3 places
- Line 500 (SUT): one formula (hardcoded MW subset, no H2)
- Line 556 (BNS in `gas_tc_pc`): a second formula including H2, using named constants
- Line 1085 (BNS viscosity in `gas_ug`): a third spelling of the same BNS formula as line 556, but indexed through `mws_lbc[0..3]`

The BNS version (556, 1085) is literally identical in intent but expressed differently. A helper `_hc_sg(sg, co2, h2s, n2, h2)` (with the `_BNS_SG_METHANE` floor optional via parameter) would unify them and remove one magic-number-indexed form. Related: `max(sg_hc, _BNS_SG_METHANE)` appears at both lines 559 and 1089 — same floor.

### M-5 — `mws` variable assigned but never read in `z_bur` (line 864)
`mws, tcs, pcs = _BNS_MWS.copy(), _BNS_TCS.copy(), _BNS_PCS.copy()` — a grep over the 80-line z_bur scope shows `tcs` and `pcs` are used (lines 875-876, 884), but `mws` has zero downstream references. Dead assignment.

### M-6 — `zee_arr, _ = convert_to_numpy(zee)` is redundant (line 1115)
In `gas_ug`, line 1034 already converts `zee` via `zee, _ = convert_to_numpy(zee)`. Line 1115 re-converts the same variable into `zee_arr` before `rhor = p / (zee_arr * R * degR * rhoc)`. `zee` is already a numpy array at that point — `zee_arr` is just an alias. Delete line 1115 and use `zee` directly.

### M-7 — Dead-equivalent duplicate computation in PMC branch of `gas_tc_pc`
Lines 481-494 compute `j, k, jt, kt` such that after line 488, `jt = j_initial + sum` and after line 490, `j = j_initial + sum` — i.e. `j` and `jt` hold identical values, as do `k` and `kt`. The code computes `tpc = kt*kt/jt` (489) and `ppc = (k*k/j)/j` (494). The result is always `ppc = tpc/j` since `k==kt` and `j==jt`. Either the original intent was to compute `jt`/`kt` using *different* coefficients (and the block should be audited against the original Piper paper), or half the arithmetic is pointless. Either way, the current form has no visible purpose and a reader cannot tell which. Action: either delete the `jt`/`kt` pair and use `j`/`k` only, or clarify via comment why two independent accumulations are preserved.

### M-8 — `gas_pvt` unpacking duplicated in both rate functions
Lines 292-296 (`gas_rate_radial`) and 376-380 (`gas_rate_linear`) are byte-identical 5-line blocks extracting sg/co2/h2s/n2/h2/zmethod/cmethod/tc/pc from `gas_pvt`. Extract to `_unpack_gas_pvt(gas_pvt)` returning a dict or tuple. Also duplicated GPvt handling in the metric branch for tc/pc (307-310 vs 387-391).

### M-9 — Nested helper functions in `gas_z` pay re-definition cost on every call
`zdak`, `z_hy`, `z_wyw`, `z_bur` are defined inside `gas_z` (lines 778, 817, 857, 870). They close over `pprs`/`tr`/`co2`/`h2s`/`n2`/`h2`/`tc`/`pc`/`is_list` so can't be trivially hoisted without a refactor, but since `gas_z` is the hottest function in the module and can be called thousands of times from `gas_dmp` / `gas_cg` / `gas_bg`, the function-object construction cost on every call is non-trivial. Either (a) lift to module-level `_zdak(pprs, tr, is_list)`, `_z_hy(...)`, etc., or (b) rely on Rust for the hot path (already the case for DAK/HY/BNS when RUST_AVAILABLE). Current structure leaves the Python fallback paying ~4× unnecessary closure setup per call.

### M-10 — `zfuncs` dispatch dict built on every `gas_z` call but used only once
Line 947: `zfuncs = {"DAK": zdak, "HY": z_hy, "WYW": z_wyw, "BNS": z_bur}`. Only one key is ever looked up (line 977/979). The dict allocation is pure overhead. Replace with a plain `if/elif` chain, or hoist the dispatch once the zfuncs become module-level (see M-9).

### M-11 — `_h2_method_override` pair-returning helper is awkward for `gas_tc_pc`
Line 472 in `gas_tc_pc`: `_, cmethod = _h2_method_override(h2, 'DAK', cmethod)`. The throwaway `_` on the zmethod side (since `gas_tc_pc` doesn't track zmethod) suggests the helper wants a companion `_h2_cmethod_override(h2, cmethod)` for single-return use, or a keyword-flag form. Minor.

### M-12 — BNS pseudocritical inner helpers `tc_ag`, `tc_gc`, `pc_fn`, `pseudo_critical` nested inside `gas_tc_pc`
Lines 527-553 define four nested functions in the `"BUR"/"BNS"` branch. All are pure, stateless, and paper-defined; they belong as module-level helpers (e.g., in a `_bns_eos.py` sub-file, see H-1). Currently re-created on every call for the ~1/3 of gas_tc_pc invocations that use BNS.

### M-13 — No `@rust_accelerated` decorator candidates identified — but many manual try/except dispatches
Scanned all 4 Rust-dispatch sites (`gas_z` 950-974, `gas_ug` 1053-1067, `gas_ponz2p` 1364-1380, `gas_dmp` 1530-1544). None can use the current `@rust_accelerated` decorator because each has (a) a method-gate (`cname in ('SUT', 'BNS')` etc.) and (b) pre-processing (unit conversion, method resolution, array shaping) that runs before the Rust call. This is a limitation of the decorator design rather than a code smell here — just noting that the current explicit-dispatch form is the correct choice for this module and should not be mechanically converted.

## Low (nitpick)

### L-1 — `darcy_gas` has no type annotations on its return's subcomponents and no validation for `l1`/`l2 > 0`
Already exported (api_ux L-1 noted its `__all__` presence). Public surface; would benefit from `validate_pe_inputs` and explicit checks for `l1 > 0`, `l2 > 0` (caller provides `r_w`/`r_ext` or `area`/`length` which the rate wrappers validate, but a direct caller is unguarded).

### L-2 — `gas_grad2sg` bisection bracket bounds `0.06958923023817742` and `3.0` are unnamed
Line 1469. The lower bound is H2's SG (≈ `MW_H2/MW_AIR` = `2.016/28.97` ≈ 0.0696). Name it `_GRAD2SG_SG_LO = MW_H2 / MW_AIR` and `_GRAD2SG_SG_HI = 3.0` with a comment on why 3.0.

### L-3 — `gas_ug` `ugz` parameter missing type annotation (line 994)
Already noted in api_ux M-12. Pure polish.

### L-4 — `_metric_to_field_pvt` does not convert `pc` back on early return, but that's fine because `tc*pc > 0` short-circuit is the only early return
No action needed; noted for future maintainers that the helper is one-way.

### L-5 — GasPVT methods have no type hints
api_ux M-11 noted this; repeating here for completeness as tech-debt polish. Trivial `-> np.ndarray` / `-> float` adds.

### L-6 — Missing test coverage for metric-unit paths of bare module functions
Only `GasPVT(... metric=True)` and `gas_hydrate(metric=True)` are exercised by test_gas.py / test_doc_examples.py. The bare functions `gas_z(metric=True)`, `gas_bg(metric=True)`, `gas_ug(metric=True)`, `gas_cg(metric=True)`, `gas_den(metric=True)`, `gas_ponz2p(metric=True)`, `gas_grad2sg(metric=True)`, `gas_dmp(metric=True)`, `gas_tc_pc(metric=True)`, `gas_rate_radial(metric=True)`, `gas_rate_linear(metric=True)`, `gas_water_content(metric=True)`, `gas_fws_sg(metric=True)` have no direct tests — any future change to any of the nine inline metric-conversion blocks could silently break without test failure. One parameterised metric round-trip test (field→metric→field consistency check) per PVT function would cover this gap.

### L-7 — `darcy_gas` has no direct unit test
`grep darcy_gas /home/mark/projects/pyResToolbox/pyrestoolbox/tests/` returns zero hits — it's only exercised transitively through `gas_rate_radial` / `gas_rate_linear`. Since it's in `__all__` (public), at least one direct test would be appropriate.

### L-8 — `gas_sg` has no test
`gas_sg(hc_mw, co2, h2s, n2, h2)` at line 1235 is public (`__all__` entry) but untested (tests/test_gas.py grep shows no usage). Trivial to add.

### L-9 — `_compute_delta_mp` has a subtle three-branch `if/elif/else` on array sizes that would read more clearly as a helper-per-shape
Lines 193-227 use `pr.size + pwf.size == 2`, `pr.size > 1`, else-branch. Readability-only — current form is correct and tested.

### L-10 — `Østergaard` cubic root solver duplicates Newton logic already in the codebase
`_required_concentration` (line 1929) is a 25-line hand-rolled Newton. The module already imports `bisect_solve` from shared_fns. For a cubic with known monotonic behaviour on (0, 100), bisection would be ~5 lines. Marginal — Newton is faster and correctness is already tested.

---

# Summary

**Counts:** 3 High, 13 Medium, 10 Low (26 items total).

**Top 3:**

1. **Split the 2195-line `gas.py` monolith** (H-1) — hydrate prediction (≈500 lines, fully self-contained) is the obvious first carve-out, then BNS/PR-EOS internals (≈130 lines) and Darcy/rate engine (≈240 lines). Oil already went through this refactor; gas has not.
2. **Extract method-resolution into a shared helper** (H-2) — the `_h2_method_override` + `validate_methods` pair is duplicated 10× across the public API, with an inconsistent BNS cross-coupling block in only 2 of the 10 sites. One `_resolve_methods()` helper eliminates the duplication and the inconsistency.
3. **Purge unused imports and harmonise metric conversion** (M-1 + H-3) — 22 dead names on the import lines, and 6 inlined duplicates of the `_metric_to_field_pvt` conversion pattern (some with over-defensive `isinstance(p, (int, float))` guards on `np.asarray * float`, which works for both scalars and arrays).
