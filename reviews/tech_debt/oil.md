# Tech Debt — oil

Scope: `pyrestoolbox/oil/` 10 sub-files (~2 900 LOC) + `tests/test_oil.py`. Tier-2 DONE items (split, named constants, `validate_pe_inputs`) not re-raised. Correctness & API/UX separately reviewed and skipped.

## High

- **Unused imports across 5 sub-files (verified by AST walk, not lint).**
  - `_compressibility.py:3` — `import numpy as np` never referenced
  - `_compressibility.py:7,8` — `SCF_PER_STB_TO_SM3_PER_SM3`, `INVBAR_TO_INVPSI` imported, never used (only PSI<->BAR direction is used)
  - `_compressibility.py:19` — `oil_deno`, `oil_viso` imported from `_correlations`, never called (dispatch goes via `oil_bo`)
  - `_correlations.py:26` — `oil_sg` imported, never used
  - `_density.py:7` — `psc` imported, never used
  - `_harmonize.py:6` — `psc` imported, never used
  - `_tables.py:3` — `import numpy as np` never referenced
  - `_tables.py:11` — `oil_sg` imported, never used
  - `_tables.py:12` — `oil_pbub` imported, never used (resolution goes through `oil_harmonize`)
  Every one is dead import. Easy delete; no behavior change.

- **`_cofb_mccain` lives in wrong sub-file.** It is a *compressibility* helper (McCain Eq 3.13 for `cofb`) defined at `_density.py:28` because `Deno_p_gt_pb` happens to use it. `_compressibility.py:18` has to reach across and import from `_density` — the exact inversion of the intended dependency direction (`_compressibility` already depends on `_correlations` which depends on `_density`). Move `_cofb_mccain` to `_compressibility.py` and have `_density.py` import it from there. If that creates a cycle, push it to `_utils.py` (it is a pure closed-form polynomial with no oil-specific helpers).

- **`__init__.py:75-76` re-exports private helpers `_cofb_mccain` and `_perrine_co_sat`.** Nothing in the package, tests, or docs imports them from the package root — grep confirms only internal usage. The top-level `__all__` already omits them. Delete both lines; keeps the surface area honest.

- **Rust `oil_bo_mccain_rust` is a dead dispatch target.** `src/oil/mod.rs:52` and `src/lib.rs:79` register it in the PyO3 module, but no Python code ever calls `_native.oil_bo_mccain_rust`. Either wire it into `oil_bo` (MCAIN branch at `_correlations.py:418-430` is the obvious target — it just calls `oil_deno` + one arithmetic line, already Rust-accelerated for SWMH), or remove the Rust function and its registration. Carrying the Rust implementation forever with no caller is net-negative maintenance.

## Medium

- **`@rust_accelerated` decorator opportunity missed for `oil_deno`.** `_density.py:187-194` uses manual `if _RUST_AVAILABLE and denomethod.name == "SWMH": try/except (ImportError, AttributeError): pass`. The `rust_accelerated` decorator was specifically introduced in Tier-3 for this pattern. The wrinkle is that the call only fires when `denomethod == SWMH`, but since that's the only enum value (see below), a plain decorator works. Worst case an enum gate can be added inside the decorator's wrapped function or the dispatch stays manual — but at minimum the imports at lines 13-15 can use the decorator's already-loaded `_rust` reference rather than a fresh conditional import.

- **`deno_method`, `co_method`, `uo_method` are single-value enums.** `classes/classes.py:81-88` — each has exactly one member (`BR`, `SWMH`, `EXPLT`). The `fn_dic = {"EXPLT": Co_explicit}` at `_compressibility.py:233` and `fn_dic = {"SWMH": ..., "PGTPB": ...}` at `_density.py:174` are dispatch tables with one correlation option. Either drop the enum + dispatch ceremony and call the implementation directly, or document "reserved for future methods" in each enum. The current state is half-built abstraction carrying per-call `validate_methods` cost.

- **Magic-number arrays in `_separator.py` not extracted.** `sg_st_gas` at `_separator.py:57-64` inlines a 5×5 Valko-McCain (2003) Eq 4-2 coefficient matrix plus a 5-term polynomial `1.219 + 0.198 Z + 0.0845 Z² + …`; `oil_rs_st` at `_separator.py:107-111` inlines a 3×3 matrix plus `3.955 + 0.83 Z - 0.024 Z² + 0.075 Z³`. Every other VALMC correlation coefficient set was extracted to `_constants.py` in Tier-3 — these two got missed. Name them `_SGST_VALMC_C`, `_SGST_POLY`, `_RSST_VALMC_C`, `_RSST_POLY`.

- **Duplicated numerical-derivative boilerplate.** `_perrine_co_sat` (`_compressibility.py:22-48`) and `Co_explicit` (`_compressibility.py:172-231`) both recompute `dp = max(0.5, p*0.001); p_hi = p+dp; p_lo = max(p-dp, psc); span = max(span, 1e-10 fallback)`. Extract to `_central_diff_bracket(p, rel=1e-3, min_dp=0.5)` returning `(p_lo, p_hi, span)`.

- **`oil_api` logic inlined at `_density.py:183`.** `api = _API_NUMER / sg_o - _API_DENOM` is exactly what `_utils.oil_api(sg_o)` does. `oil_api` is already re-exported publicly. Same import line already pulls `oil_sg`; just add `oil_api`.

- **29 `if metric:` blocks across the 10 sub-files.** Each public function repeats an ~8-line conversion stanza (`p *= BAR_TO_PSI`, `degf = degc_to_degf(degf)`, conditional rsb/pb/p_uo multiplications). `_rate.py` additionally wraps those in a useless `if not isinstance(pr, (int, float))` guard (`np.asarray(scalar)` already handles scalars — lines 51-53, 140-142 are dead defense). A `_metric_to_field(**kwargs)` helper in `_utils.py` that converts a dict of labelled values would cut ~200 lines and centralize the unit contract.

- **OilPVT `rs()`/`bo()`/`density()`/`viscosity()` share a 10-line recursion+metric+rs-default prelude.** `_pvt_class.py:100-108`, `110-129`, `131-153`, `155-172` — four methods with the same `if self._is_array(p): recurse; if self.metric: p_field/degf_field; if rs is None: _rs_field else convert` opener. Pull it out into a `_prepare_point(p, degf, rs)` helper returning `(p_field, degf_field, rs_field, is_array)`.

- **Test coverage gaps by sub-file.**
  - `_utils.py`: `oil_ja_sg()` has zero direct tests (only exercised transitively through `oil_twu_props`); `get_real_part` has zero tests and the function is trivial enough it should probably inline at the 2 call sites in `_correlations.py` (lines 306, 339) and be deleted.
  - `_separator.py`: `sg_st_gas()` and `sgg_wt_avg()` have zero tests. These are publicly re-exported in `__init__.py:54-57` and in docstrings.
  - `_tables.py`: `_resolve_pb_rsb`, `_build_bot_tables`, `_format_bot_results` are only tested transitively via `test_make_bot_og`. The 100-iteration PVTO extension loop at lines 45-58 has no direct coverage. If extraction was worth the split, unit tests should follow.
  - `_harmonize.py`: the `RuntimeError("Could not solve …")` path at line 113 has no test (API review noted the error-case gap separately; this is the harmonize-path variant).
  - `_compressibility.py`: `_perrine_co_sat` has no direct tests; only exercised through `oil_co(co_sat=True)`.

## Low

- **`check_sgs` unused in some call sites.** `oil_rs` (`_correlations.py:292`) calls `check_sgs(sg_g=0, sg_sp=sg_sp)` — with `sg_g=0` hardcoded. The function always returns `(sg_g, sg_sp)` with `sg_g` imputed from `sg_sp`, yet `sg_g` is never used downstream in `oil_rs` (Velarde only needs `sg_sp`, Valko-McCain only needs `sg_sp`, Standing uses `sg_g`). When `sg_g` is not needed the call is dead work; when it is (Standing path), the user never supplied it. Minor, but the dead branch should be noted.

- **`get_real_part` is one-line pass-through.** `_utils.py:10-14` — 5 lines to return `value.real if complex else value`. `np.real(value)` does the same and is already imported via `np`. Replace with `np.real`, delete helper, drop from imports.

- **`_utils.py` bucket: `oil_twu_props` is not a "utility".** It is a 160-line correlation-heavy function with nested helpers for paraffin properties, boiling point iteration, and a Twu(1984) recipe. `_utils.py` should hold the ~5 trivial helpers (SG↔API, `check_sgs`, `get_real_part`). `oil_twu_props` belongs in its own `_twu.py` or merged into `_correlations.py`. Currently 160 of `_utils.py`'s 226 lines are one function.

- **Duplicated imports inside `oil_rs_bub`.** `_correlations.py:209,239` do `import warnings` inside function bodies even though `warnings` is imported at module level (line 3). Dead reimport.

- **Type hint gaps.**
  - `oil_co` (`_compressibility.py:51-70`): no return annotation. The function returns `float` OR `list[float]` depending on `co_sat`, which actually deserves documentation via `Union[float, list[float]]`.
  - `OilPVT.__init__`, `rs`, `bo`, `density`, `viscosity`: zero type annotations. Inconsistent with the rest of the module (module-level functions carry full annotations).
  - `_rate.py:oil_rate_radial/linear`: `oil_pvt` parameter typed as `None` rather than `Optional[OilPVT]`.
  - `_separator.py:sgg_wt_avg` has no return-type annotation.

- **`_separator.py:sg_evolved_gas` / `sg_st_gas` are inconsistent with module convention.** No `metric` parameter, no `validate_pe_inputs` call. Every other public oil function accepts metric input; these two are field-only. Either document the restriction in the docstring or add the conversion for consistency.

- **`oil_harmonize_pb_rsb` + `OilPVT.from_harmonize` carry no `DeprecationWarning`.** API-UX review already flagged these, noting here only to confirm: once warnings are added, the deprecated paths can be dropped on next major, at which point the `_pvt_class.py:72-84` classmethod shrinks to 0 lines.

- **`pbmethod = 'VALMC'` shadowed inside `oil_rs_bub`.** `_correlations.py:172-175` sets `pbmethod = 'VALMC'` as a local string, then `validate_methods(["pbmethod", "rsmethod"], [pbmethod, rsmethod])` overwrites it with the enum. Only the VALMC branch (`rsbub_valko_mccain`) actually uses `pbmethod`, and it's always VALMC there by construction. Either drop the `pbmethod` local entirely or pass it as an explicit constant to the inner function. Mild confusion at a glance.

- **Inner-function default-arg closure trick in `Co_explicit`.** `_compressibility.py:172-185` declares defaults `p=p, api=api, …` to "freeze" values. Works, but is a Python idiom that trips reviewers. A named helper at module scope taking explicit args is clearer.

---

## Summary

Counts: 4 High + 8 Medium + 9 Low tech-debt items. Top 3:

1. **Confirmed dead imports in 5 of 10 sub-files** (numpy unused in `_compressibility.py`/`_tables.py`, `oil_sg`/`oil_pbub`/`oil_viso`/`oil_deno`/`psc`/`SCF_PER_STB_TO_SM3_PER_SM3`/`INVBAR_TO_INVPSI` dead across multiple files) — trivial to delete.
2. **`_cofb_mccain` lives in `_density.py` but belongs in `_compressibility.py`**, forcing a reverse import; the `__init__.py` re-exports of private helpers `_cofb_mccain` and `_perrine_co_sat` (lines 75-76) are dead.
3. **Dead Rust dispatch path** (`oil_bo_mccain_rust` registered but never called), plus single-value "enums" (`deno_method`, `co_method`, `uo_method`) driving speculative `fn_dic` dispatch that never materialized.
