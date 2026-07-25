# API/UX — matbal

## Critical

### 1. No aquifer model plug-in — `We` must be precomputed by user
`gas_matbal` accepts `We` as a pre-computed array. `oil_matbal` has no `We` input at all — water influx cannot be incorporated into the oil MB. The standard workflow (couple matbal to a Fetkovich/Carter-Tracy aquifer and iterate) is unsupported. `simtools.influence_tables` exists for VEH tables but is not linked into matbal. RST claims the module supports "Havlena-Odeh aquifer influx correction", but there is no aquifer model shipped — the user must bring their own `We` time series.

**Impact**: Primary real-world use case (water-drive reservoirs) requires off-library code. No closure on a canonical petroleum-engineering workflow.

### 2. `oil_matbal` has no `We` / aquifer term whatsoever
Havlena-Odeh oil MB formally is `F = N(Eo + mEg + (1+m)Efw) + We`. `oil_matbal` omits `We` entirely — there is no parameter to pass cumulative aquifer influx into the oil calculation. Users with water-drive oil reservoirs cannot use this function correctly. This is a missing API surface, not just a missing aquifer model.

## High

### 3. No plotting hooks / no matplotlib figure returned
Cole plot (`cole_F_over_Et`), P/Z vs Gp, Havlena-Odeh straight-line plot, drive-index stacked plot — none of these canonical MB diagnostic plots are produced. The RST shows no plotting example. `dca` and `nodal` modules have notebook examples with plots; `matbal` does not appear in `examples.ipynb` patterns for visualisation. Users must hand-roll every plot from the dataclass arrays. At minimum a `plot()` method on the result dataclasses, or a helper in `matbal` that returns a `matplotlib.figure.Figure`, would match typical MB tool UX (MBAL, IHS Harmony).

### 4. No DataFrame input path
Real matbal inputs typically come from a well history DataFrame (date, p, Gp/Np, Wp, etc.). Both functions take `array-like` only and require the user to split columns manually. A pandas DataFrame + column-name-mapping dict convention would cut boilerplate considerably and align with how engineers actually hold this data. Not blocking but painful.

### 5. `regress` API — dict-of-tuples is awkward, error messages inconsistent
`regress={'m': (0, 2), 'cf': (1e-7, 50e-6)}` is workable but:
- No way to fix initial guess independently of bounds (the function picks midpoint if current value is outside bounds — silently).
- Allowed keys hard-coded to `{'m','cf','cw','sw_i'}`. Error `Invalid regress key 'X'. Allowed: [...]` is good; but `regress['X'] must be a (lower, upper) tuple` uses backticks vs single quotes inconsistently with other messages.
- No way to regress OOIP (N) directly — the function always optimises CV of per-step N estimates; some workflows prefer to regress N plus aquifer params together.
- No optimiser choice exposed (always L-BFGS-B); no convergence diagnostics returned (no `success`, `nit`, `fun` from scipy result).

### 6. `pvt_table` validation message too terse
`"Survey pressures must fall within pvt_table pressure range"` does not report the offending pressure or the table range. Engineers debugging a 20-survey sequence will have to inspect arrays by hand. Should report `"Survey pressure {p_bad} outside pvt_table range [{pmin}, {pmax}]"`.

### 7. Gas `pvt_table` schema inconsistent with oil `pvt_table`
Gas: `{'p', 'Z'}` OR `{'p', 'Bg'}` (either). Oil: `{'p', 'Rs', 'Bo', 'Bg'}` (all four required, no options). No documentation of why the gas case accepts two schemas but oil accepts only one. Users migrating between functions will stumble.

## Medium

### 8. `Gp` units are silently overloaded
Docstring says Gp units are "user-defined" so OGIP comes back in the same units — but then says "When Wp or We are provided, Gp should be in scf (or sm3 if metric) for dimensional consistency with Bg". This is a load-bearing caveat buried in the Gp parameter description. Mixing arbitrary Gp units with real-volume Wp/We terms silently produces wrong F. Should be enforced with a warning, or surfaced in a dedicated `.. warning::` block in RST, not a single sentence in a bullet.

### 9. Enum/string method dispatch — no validation at entry
`zmethod='DAK'`, `cmethod='PMC'`, `rsmethod='VELAR'`, `bomethod='MCAIN'` accept strings but there's no call to the `validate` module (as required by project conventions for new code). A typo like `zmethod='DACK'` will fail deep in `gas.gas_z` with a less-helpful error.

### 10. Drive-index dict nesting — awkward access pattern
`r.drive_indices['DDI']` requires two lookups. Either three attributes on the dataclass (`r.ddi`, `r.sdi`, `r.cdi`) or an `Nx3` numpy array with a column labels property would be friendlier. Current form mixes flat attribute access (F, Eo, Eg, Efw) with nested dict access (drive_indices), which is inconsistent.

### 11. Missing `validate_pe_inputs` at entry
Per Tier-1 "New Code Requirements", public entry points must call `validate_pe_inputs` after metric conversion. Neither `gas_matbal` nor `oil_matbal` does. Inputs like `sg=0` or `co2=1.5` are not caught here — they'd fail in `gas_z` with a harder-to-trace error.

### 12. `OilMatbalResult.regressed` is `dict | None`, no indication of convergence
When `regress` is used, user gets back `{'m': 0.42, 'cf': 4.7e-6}` with no info on whether L-BFGS-B converged, how many iterations, final objective value, or whether optimum hit a bound. Hitting a bound is a strong diagnostic signal that the regression is meaningless — currently invisible.

### 13. No warning when P/Z R-squared is poor
`r_squared` is returned but not checked. A value < 0.9 usually indicates aquifer influx or bad data. No user-facing warning; user must inspect manually.

## Low

### 14. Attribute naming inconsistency — `gp` (lowercase) vs `F`, `Et` (PascalCase)
`GasMatbalResult` has `.gp` (lowercase) alongside `.F`, `.Et`, `.bg`. Input is `Gp` (PascalCase). Mixed conventions within one dataclass. Minor but jarring.

### 15. `method` field in `GasMatbalResult` is a string, not an Enum
`method='pz'` or `'havlena_odeh'`. Elsewhere the library uses Enum dispatch for method selection. A string return is acceptable here (it's output not input) but mildly inconsistent.

### 16. No `__repr__` customisation on result dataclasses
Printing a `GasMatbalResult` dumps all arrays including full `pz`, `gp`, `bg`, `F`, `Et`, `cole_F_over_Et` numpy arrays. A concise repr showing ogip, r_squared, n_points, method would improve notebook/REPL UX.

### 17. Docstring/RST drift — `Bg` units in RST
RST says `bg` units are `rcf/scf | rm3/sm3`. Docstring in dataclass says only `rcf/scf`. Both are technically right (rcf/scf == rm3/sm3 numerically) but inconsistent wording.

### 18. `GasMatbalResult.method` value `'havlena_odeh'` uses snake_case; dataclass docstring also snake_case; but `OilMatbalResult` has no analogous method tag
Asymmetric. Oil MB has only one internal path; harmless, but if aquifer support is added (Issue 2), this will need parity.

### 19. No `to_dict()` / `to_dataframe()` on result classes
Common ask for downstream reporting. One-liner methods would help.

### 20. Compressibility guidance block is in RST only, not in docstring
The `.. note::` about `cw_usat` vs `cw_sat` is in `matbal.rst` but not in the `oil_matbal` docstring. Agentic users reading docstrings at runtime miss it.

---

## Summary

**Counts**: 2 Critical, 5 High, 6 Medium, 7 Low.

**Top 3**:
1. **No aquifer model plug-in and oil_matbal has no `We` at all** — the flagship petroleum-engineering workflow (water-drive MB with Fetkovich/Carter-Tracy) is absent. `simtools.influence_tables` is orphaned.
2. **No plotting hooks** — Cole plot, P/Z plot, Havlena-Odeh straight line, drive-index stack are not produced. Users hand-roll every plot. A `.plot()` method or helper in the module would match tool-category UX.
3. **No DataFrame input and weak error messages** — real inputs come from DataFrames; array-of-columns boilerplate plus terse `pvt_table` out-of-range errors without offending value create friction.
