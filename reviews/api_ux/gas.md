# API/UX — gas

Scope: public API of `pyrestoolbox/gas/gas.py` and alignment with `pyrestoolbox/docs/gas.rst`. Correctness/numerical findings already covered in `reviews/correctness/gas.md` are NOT repeated here.

---

## Critical (breaks user workflow / wrong result from confusing UX)

None found.

---

## High (inconsistency, misleading default, silent type coercion)

### H-1 — `zmethod='BUR'` is undocumented in code docstrings yet is the primary example in the RST
Every in-code docstring (`gas_z`, `gas_ug`, `gas_bg`, `gas_cg`, `gas_den`, `gas_ponz2p`, `gas_grad2sg`, `gas_dmp`, `gas_rate_radial`, `gas_rate_linear`) lists the accepted string values as `'DAK' | 'HY' | 'WYW' | 'BNS'`. But the RST (`docs/gas.rst` lines 25, 32, 71, 79) uses `'BUR'` prominently, including in the canonical pure-CO2 example. Because `z_method.BNS = z_method.BUR = 3`, both *work*, but a user reading only the docstring sees `BNS`, while a user reading the RST sees `BUR`. This is discoverability-debt: either decide `BNS` is canonical and purge `BUR` from the RST, or list both in docstrings. The Tier 1 CLAUDE.md note "BNS canonical, BUR alias" suggests `BNS` is canonical; the RST should be updated to match (primary examples use `BNS`, with a parenthetical "`BUR` is accepted as a legacy alias").

### H-2 — `c_method` enum has both `BUR` and `BNS` as distinct members, but they dispatch identically
`classes.py:64–65` defines `c_method.BUR = 2`, `c_method.BNS = 3` as two *different* enum values (unlike `z_method` where BNS and BUR share value 3). Yet `gas_tc_pc` line 526 handles both via the same branch `if cmethod.name in ("BUR", "BNS")`. This is invisible to string users (validate_methods maps `'BUR' → 'BNS'` via `legacy_map`, so string dispatch always returns `c_method.BNS`), but a user who imports the enum and passes `c_method.BUR` directly gets enum value 2, which then works correctly only because of the `.name in (...)` check. Either (a) collapse to a single alias like `z_method` does, or (b) document the distinction. Current state is confusingly asymmetric across the two enums.

### H-3 — `gas_ug` `zee` parameter silently recomputes on mismatched length without warning
`gas.py:1032–1048`: when a user provides `zee` as a scalar but `p` as an array (or vice-versa) `check_2_inputs` returns False and `zee` is silently replaced by an internal `gas_z(...)` call. The user believes they supplied `zee` but it was discarded. At minimum emit a `warnings.warn` when a non-default `zee` is dropped; ideally raise `ValueError(f"zee length {len(zee)} does not match p length {len(p)}")`. Also: the `zee_provided` sentinel logic `not (isinstance(zee, (int, float, bool)) and not zee)` is opaque — a user passing `zee=0.0` (legitimate near-zero guess) would be considered "not provided", which is probably the intent but worth a docstring note.

### H-4 — `gas_rate_radial` / `gas_rate_linear` skip `sg` validation when `gas_pvt` supplied, but require it otherwise — silent divergence on invalid `gas_pvt`
`_prepare_gas_rate_inputs` (line 231) calls `validate_pe_inputs(..., sg=sg, ...)`. When `gas_pvt` is provided, `_prepare_gas_rate_inputs` is bypassed (line 319/398), so a malformed `GasPVT` object with e.g. `sg=0` or unset composition is never revalidated. `GasPVT.__init__` does not call `validate_pe_inputs` at all — users can construct `GasPVT(sg=-0.5)` without error. Add `validate_pe_inputs` in `GasPVT.__init__`.

### H-5 — `gas_grad2sg` docstring says "fail if gas SG is below 0.55, or greater than 3.0" but bisection upper bound is 3.0 while RST claims 1.75
Docstring (`gas.py:1417`) says "below 0.55 or greater than 3.0". RST (`docs/gas.rst:733`) says "below 0.55, or greater than 1.75". Actual code: bracket is `[0.06958923023817742, 3.0]` (line 1469). Three inconsistent sources. The lower bound (0.0696) is also bizarre (hydrogen has SG ≈ 0.0696; the comment should explain this is the H₂ floor). Align all three to reflect actual bisection bracket. The lower bound of ~0.07 is not "0.55" — the docstring overstates the floor.

### H-6 — `GasPVT` silently accepts `metric=True` without validating that `tc`/`pc` are pre-converted
`GasPVT.__init__` never accepts `tc`/`pc` user overrides (unlike `gas_tc_pc()`), so all critical properties come from `gas_tc_pc(sg, co2, h2s, n2, h2, cmethod.name)` called at line 1638 in oilfield units. This is fine, but undocumented: a user reading the constructor signature sees `metric=False` and might expect to be able to pass their own metric-unit `tc`/`pc`. Either add explicit `tc=0, pc=0` parameters (with metric-aware conversion) to mirror `gas_tc_pc`, or document that user overrides of Tc/Pc are not supported by `GasPVT` and to use bare `gas_z` instead.

---

## Medium (naming, docstring gaps)

### M-1 — Docstring typo: `h2` in `gas_bg` says "Nitrogen"
`gas.py:1221`: `h2: Molar fraction of Nitrogen. Defaults to zero if undefined` — should be Hydrogen. (This is already flagged in correctness/gas.md as "Low" but it is a pure docstring issue; listing here too for API/UX completeness.)

### M-2 — `gas_dmp` parameter order is `(p1, p2, degf, sg, ...)` while every other gas function uses `(p, sg, degf, ...)`
All other PVT functions (`gas_z`, `gas_ug`, `gas_bg`, `gas_cg`, `gas_den`, `gas_ponz2p`, `gas_grad2sg`) use the ordering `(p, sg, degf, ...)`. `gas_dmp` uses `(p1, p2, degf, sg, ...)` — `sg` and `degf` are swapped relative to the family convention. Not a bug, but keyword-only callers are shielded; positional-args callers hit a surprise. Cannot be changed without breaking compat; note in the docstring.

### M-3 — Inconsistent argument order: `gas_cg` vs other `gas_*` functions
`gas_cg` signature (line 1124): `(p, sg, degf, co2, h2s, n2, h2, tc, pc, zmethod, cmethod, metric)` — impurity mole fractions come BEFORE `zmethod`/`cmethod`. Every other function places method args immediately after the required inputs: `(p, sg, degf, zmethod, cmethod, co2, ...)`. This means a keyword-free positional call like `gas_cg(2000, 0.75, 180, 'DAK')` errors out in surprising ways (`'DAK'` lands on `co2` → string-as-mole-fraction). Code in repo uses keyword args, but downstream positional users will be burned. Low likelihood of it mattering in practice but worth a note.

### M-4 — `gas_sg` has no `validate_pe_inputs` call; silently accepts nonsense inputs
`gas_sg(hc_mw, co2, h2s, n2, h2)` (line 1235) accepts any values and returns a float. Negative `hc_mw` or mole fractions summing to >1 return garbage SG. Per CLAUDE.md "New Code Requirements", all public functions should validate. Add `validate_pe_inputs(co2=co2, h2s=h2s, n2=n2, h2=h2)` and a check that `hc_mw > 0`.

### M-5 — `gas_fws_sg` has no `validate_pe_inputs` and no docstring for `api_st` range
`gas_fws_sg(sg_g, cgr, api_st, metric)` (line 1587) silently accepts negative API, negative CGR, negative SG. Add input validation. Also docstring order mismatches signature: signature is `(sg_g, cgr, api_st)`, docstring lists `sg_g`, `api_st`, `cgr` (lines 1593–1595) — a reader cross-referencing the two gets tripped up.

### M-6 — `gas_water_content` lacks upper salinity bound / warning
`gas_water_content(p, degf, salinity=0, metric=False)`: the Danesh salinity correction is `(1 − 0.00492·s − 0.00017672·s²)`. For large salinity (say 30 wt%), this becomes negative, producing negative water content. Add a validity-range warning (Danesh correlation is calibrated to ~25 wt% NaCl) or raise for `salinity >= 20`.

### M-7 — `gas_grad2sg` param table docstring describes `p` as "float, list or np.array" in RST (line 746), but code only accepts scalar
`gas.py:1402–1414`: `p: float`. RST says `float, list or np.array`. Bisection `bisect_solve` is scalar. Correct the RST to `float`.

### M-8 — `gas_hydrate` returns `hfp = NaN` silently when target T is outside correlation range
`_hydrate_formation_press` (line 1889) returns `float('nan')` if T is outside `[t_lo, t_hi]`. This is then quietly passed through to `HydrateResult.hfp`. A user plotting `hfp` sees a missing value without knowing why. Either log a `warnings.warn` or explicitly document in the `HydrateResult.hfp` docstring that NaN means "out of correlation range" (currently only mentioned cryptically as "or NaN if no inhibitor" in a different field's doc).

### M-9 — `metric` not propagated into reported error messages
`_ponz2p_err_msg(target)` (line 1354) always prints "P/Z={target:.4f}" using the *converted* (psia) value. If a metric-units user gets an error, the reported number is in psia, not barsa — surprising. Either convert back on error, or append "(psia, internal units)".

### M-10 — `gas_dmp` `p1 == p2` short-circuits to return `0` without type (float vs array)
`gas.py:1523–1524`: returns bare `0` (int). Other branches return float. Minor typing wart.

### M-11 — Missing return-type annotations on `GasPVT` methods
`GasPVT.z`, `.viscosity`, `.density`, `.bg` have no type hints (lines 1647–1676). Other public functions have `-> np.ndarray` or `-> float`. Trivial polish.

### M-12 — `ugz: bool = False` in `gas_ug` signature has no type annotation
Line 994: `ugz = False,` — no annotation. Line 1022 docstring type is "boolean". Trivial.

### M-13 — Docstring for `gas_cg` has copy-paste artefact `pwf: BHFP (psia)`
Line 1143: `pwf: BHFP (psia)` — `gas_cg` takes no `pwf` parameter.

### M-14 — `gas_rate_radial`/`gas_rate_linear` RST tables drop `h2` from the parameter list
RST lines 980–988 (radial) and lines 1138–1146 (linear) describe `n2`, `co2`, `h2s` but do NOT include `h2` or `tc`/`pc` (for linear). The code accepts all of these. RST tables are incomplete.

### M-15 — `HydrateResult.inhibited_hft` and `required_inhibitor_wt_pct` semantics overlap confusingly
`inhibited_hft` is HFT minus the *applied* dose depression; `required_inhibitor_wt_pct` is the dose *needed* to reach operating temperature. When `inhibitor_wt_pct=0` but `inhibitor_type` is set, `inhibited_hft == hft` (line 2124) and `required_inhibitor_wt_pct > 0` — user may think dose is sufficient because `inhibited_hft > degf`, whereas `required` says otherwise. Docstring already notes this (line 1795–1801); consider surfacing it more prominently as a short "comparing fields" note in the function-level docstring.

### M-16 — `gas_hydrate` accepts `inhibitor_wt_pct` up to 100 (exclusive) with no validation against `_MAX_WT_PCT` for applied dose
Line 2044 checks `inhibitor_wt_pct < 0 or inhibitor_wt_pct >= 100` but does not warn if the *applied* dose exceeds the physical maximum for the selected inhibitor type (e.g., MEOH max 25%). A user passing `inhibitor_type='MEOH', inhibitor_wt_pct=80` gets an `inhibited_hft` computed from the Ostergaard cubic extrapolated far outside its range with no warning. Add `warnings.warn` when `inhibitor_wt_pct > _MAX_WT_PCT[inh]`.

---

## Low (polish)

### L-1 — `darcy_gas` is re-exported in `__all__` but is an internal helper
`__all__` line 54 lists `darcy_gas`. It is documented in neither the RST nor the module header. If intentional (user can call it directly), add RST coverage; otherwise prefix with `_` and remove from `__all__`.

### L-2 — `z_method` docstring in RST uses legacy phrase "Slowest, Most Accurate" / "Fastest, Least Accurate"
`docs/gas.rst:22–24`: subjective performance labels without benchmarks. Consider removing or backing with the current benchmark numbers (CLAUDE.md review notes a 35× disparity).

### L-3 — `gas_tc_pc` docstring line 443 says "deg R" return when `metric=False` but no unit for metric case in the docstring (the RST covers it)
Minor doc consistency.

### L-4 — Inhibitor string accepted values in `gas_hydrate` not cross-validated with module re-exports
`inhibitor` enum is re-exported (`__all__` line 57) but `gas_hydrate` accepts strings like `'MEOH'` directly via `validate_methods`. Users cannot pass `inhibitor.MEOH` (enum) — `validate_methods` only branches on `isinstance(variables[m], str)`. Enum inputs pass through unchanged, which happens to work because `inh` is used as a dict key. Worth documenting: "Accepts `inhibitor_type='MEOH'` or `inhibitor_type=inhibitor.MEOH`."

### L-5 — `Tuple` imported from `typing` but only used for one return annotation
Line 64: `from typing import Tuple`. Could use lowercase `tuple` with Python 3.9+. The project targets cp38+ so leave as-is.

### L-6 — `_metric_to_field_pvt` applied inconsistently
`gas_z`, `gas_bg`, `gas_cg`, `gas_den` call the helper. `gas_ponz2p`, `gas_grad2sg`, `gas_dmp`, `gas_rate_*` duplicate the conversion inline. Already noted in tech-debt contexts — leaving as Low here as it's a consistency nit for readers comparing function bodies.

### L-7 — RST example for `gas_rate_linear` GasPVT overload shows positional order: `(k, area, length, pr, pwf, degf, ...)`
But actual signature is `(k, pr, pwf, area, length, degf, ...)` (line 329). The RST example (line 1172, 1175, 1183) uses keyword args so it works, but new users copy-pasting may misread the order.

### L-8 — `HydrateResult` dataclass has no `__repr__` customization
A user printing a HydrateResult gets all 17 fields on one line — difficult to read. Adding a custom `__str__` summarizing key fields (hft, subcooling, required wt%) would dramatically improve interactive UX. Low priority.

### L-9 — No re-export of `HydrateResult`/`GasPVT` at the top-level `pyrestoolbox` package level is documented
The gas `__init__.py` uses `from .gas import *` which picks up `__all__`, so `from pyrestoolbox.gas import GasPVT` works. A quick doc note confirming both `from pyrestoolbox import gas; gas.GasPVT(...)` and the submodule import paths would help LLM agents.
