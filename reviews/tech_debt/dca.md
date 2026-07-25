# Tech Debt — dca

Module is a single 968-line file with eight `_fit_*` helpers, three public `fit_*` entry points, two forecast functions, and three dataclasses. Unit-agnostic (so `validate_pe_inputs` correctly does not apply). Rust acceleration only exists for the two hyperbolic fitters. No separate sub-files.

## High

### H1. Six near-identical fitter bodies — R² / residual / DeclineResult boilerplate duplicated
`_fit_exponential`, `_fit_harmonic`, `_fit_hyperbolic`, `_fit_exponential_cum`, `_fit_harmonic_cum`, `_fit_hyperbolic_cum` (lines 337-565) all repeat the same 8-line block:

```python
ss_res = np.sum((q - q_pred) ** 2)
ss_tot = np.sum((q - np.mean(q)) ** 2)
r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
return DeclineResult(method=..., qi=qi, di=di, b=..., r_squared=r2, residuals=q - q_pred)
```

A single helper `_build_decline_result(method, qi, di, b, q_obs, q_pred)` would eliminate ~50 lines. The two hyperbolic grid-search bodies (lines 376-425 and 502-565) are structurally identical — same `best_r2`, same `np.arange(0.05, 0.96, 0.01)` loop, same R² computation — only the linearization transform differs. Could collapse to one `_hyperbolic_grid_search(t_or_Np, q, transform_fn, recover_params_fn)` taking two callables, or keep separate but factor out the loop scaffold.

Similarly `fit_decline` (lines 672-737) and `fit_decline_cum` (lines 568-669) duplicate:
- windowing block (`mask`, slicing, `len == 0` guard, shift-to-origin — lines 601-613 vs 696-706)
- length / positivity validation (lines 615-620 vs 708-713)
- dispatch table + `'best'` loop + `max(results, key=r_squared)` (lines 625-645 vs 715-737)

About 40 lines of duplication across the two public fitters. A shared `_run_dispatch(fitters, x, q, method)` would also give uniform error messages (see #13 in API review).

### H2. Arps forecast vs EUR share analytic inverse — duplicated in `eur()` and nowhere else
`eur()` (lines 300-334) re-derives the time-at-qmin inverse for each b branch. Same three branches appear in `arps_rate`, `arps_cum`. Not large, but a private `_t_at_rate(qi, di, b, q_target)` helper would DRY the branch logic and make it reusable by `forecast()` for more accurate `q_min` truncation (forecast currently uses discrete `np.where(mask)[0][-1]` which is step-size dependent).

## Medium

### M1. Magic numbers unextracted — no named constants at all
The module has **zero** module-level named constants despite several correlation-like magic values:

| Line | Value | Meaning |
|------|-------|---------|
| 293 | `0.001`, `500`, `10` | Duong trap integration lower bound, min grid points, points-per-unit-t |
| 398, 527 | `0.05, 0.96, 0.01` | Hyperbolic b-grid start/stop/step |
| 442 | `[0, 0.01, 1.001], [q_f[0]*5, 10.0, 3.0]`, `maxfev=5000` | Duong `curve_fit` bounds + iters |
| 443 | `p0=[q_f[0], 1.0, 1.2]` | Duong initial guess |
| 514, 553 | `1e-10` | Clamp floor for `inner` in cum hyperbolic |
| 536 | `1e-30` | Near-zero divisor guard |
| 796-800 | `1.5`, `0.01`, `10.0`, `1e-8`, `1e-3`, `1e6`, `maxfev=5000` | Logistic ratio curve_fit |
| 925 | `dt / 2` | `np.arange` upper-bound fudge |

Per CLAUDE.md's "Named Constants" rule, these should be module-level `_ALL_CAPS` with paper citations (Arps 1945 for b-bounds, Duong 2011 for Duong bounds/initial). Current state violates the "all new code" requirement — the rule existed at creation time per the Tier-3 review history. 17 magic numbers is not trivial.

### M2. Hyperbolic grid resolution hardcoded — no way to tune speed/accuracy
`np.arange(0.05, 0.96, 0.01)` = 91 trials × full RANSAC regression each. For the two grid searches combined, this is the hot path when `method='best'`. No `b_min`, `b_max`, `n_steps` knobs, no adaptive refinement around the best point (a two-stage coarse-then-fine grid would give better b resolution than 0.01 with fewer trials). Also excludes the physically valid tails b ∈ (0, 0.05) and (0.95, 1) — see correctness review #4. The tails could be addressed simultaneously with the grid-resolution tech-debt fix.

### M3. Sentinel `0.0` pattern in `DeclineResult` — typed but semantically ambiguous
`DeclineResult` has `di=0.0, b=0.0, a=0.0, m=0.0` defaults. For Arps results, `a` and `m` are meaningless zeros; for Duong, `di` and `b` are meaningless zeros. User code reading `result.di` on a Duong result sees `0.0` (a plausible physical value) rather than "not applicable". Fixes in order of effort:
1. Use `Optional[float] = None` for method-specific fields — user sees explicit `None`.
2. Two dataclasses `ArpsResult(method, qi, di, b, ...)` and `DuongResult(method, qi, a, m, ...)` sharing a `_DeclineResultBase` with method/r_squared/residuals. `fit_decline` returns the Union. More invasive (breaks anyone constructing `DeclineResult` directly — which `forecast()` tests do).
3. `params: dict` attribute replacing positional fields. Most flexible, least type-safe.

Same problem in `RatioResult.c=0.0` (only meaningful for logistic). Prefer (1) as minimum fix, matches the API-UX reviewer's recommendation.

### M4. Missing type hints on every function signature
Public and private signatures have **no** type annotations. Example: `def arps_rate(qi, di, b, t):`, `def fit_decline(t, q, method='best', t_start=None, t_end=None):`. Docstrings describe types, but mypy / IDE autocomplete get nothing. Given the module exposes structured dataclasses (already typed) and mixes scalar/array input (ArrayLike return), adding `from numpy.typing import ArrayLike` and annotating would cost ~30 one-line edits. The rest of the package has partial type hints — DCA appears to be among the weakest.

### M5. `__all__` and `__init__.py` convention inconsistent with oil/
`dca/__init__.py` just re-exports via `from .dca import *`. Every other multi-file module (oil, simtools, etc.) has explicit re-exports. Minor, but if the file is ever split per the Oil-style refactor, the wildcard import hides the public API.

## Low

### L1. `_fit_duong` uses unbound upper bound for `qi` (`q_f[0] * 5`)
Bounds `[0, 0.01, 1.001], [q_f[0] * 5, 10.0, 3.0]` — the `q_f[0] * 5` upper bound on qi is arbitrary; for noisy data where q_f[0] is an outlier, this can silently cap the fit. Same pattern in logistic (`Rmax_guess * 5`). Magic multiplier, no citation.

### L2. `_fit_duong` variable `q_pred_full`, `q_pred_valid` shadow-redundancy
Lines 446-448 create `q_pred_full = np.full_like(q, np.nan)`, fill masked entries, then re-slice `q_pred_valid = q_pred_full[mask]`. The full-sized array is never returned — immediately sliced back to the masked view. Could be one line: `q_pred_valid = duong_func(t_f, qi, a, m)`. Minor dead computation.

### L3. `fit_decline` / `fit_decline_cum` unused `name` loop variable
Lines 633-636 and 723-727: `for name, fitter in fitters.items():` — `name` is never used inside the loop. `fitters.values()` would be clearer, or keep `name` and add it to an error-message candidate list (links to API-UX issue #8, "no tie-break reporting").

### L4. `duong_cum` per-point integration instead of shared fine grid
Lines 291-295 — loops over each `ti`, creating a new `np.linspace(0.001, ti, ...)` and re-evaluating `duong_func` from scratch. For an array of monotone t, a single fine grid + cumulative `np.trapezoid` would be O(N) instead of O(N·500). Performance tech debt, low impact for small t arrays.

### L5. RANSAC tolerance / iterations not exposed
`ransac_linreg` is called 11 times across the module with default kwargs. No way for a user to tighten outlier sensitivity or request deterministic (non-RANSAC) regression. Adds an unavoidable layer of stochasticity to every `fit_*` call — test failures have been worked around with wide tolerance ranges (e.g., `0.2 < result.b < 0.9` in `test_fit_hyperbolic_outliers`). Parameter to opt out for reproducibility / unit-test stability would be nice.

### L6. Test coverage gaps
`test_dca.py` (63 tests) is substantial but misses:
- `fit_ratio_cum` roundtrip (no test asserts `ratio_forecast(fit_ratio(...), x) ≈ ratio`).
- `duong_cum` at small `t` (<1) — the negative-value bug in correctness review #1 would have been caught.
- `forecast` with `q_min > qi` (edge case in correctness review #14).
- `forecast` with `dt=0` or negative (correctness review #2).
- `uptime` out-of-range values (correctness review #3).
- `fit_decline_cum` with non-monotone `Np` (correctness review #6).
- No Rust-vs-Python parity test for `fit_hyperbolic_rust` / `fit_hyperbolic_cum_rust` — the try/except silent fallback (lines 391, 520) means a regression could silently slip to the slower Python path with no test catching the difference.

### L7. `convert_to_numpy` / `process_output` inconsistency
`arps_rate`, `arps_cum`, `duong_rate`, `duong_cum`, `ratio_forecast` accept scalar/list/array. `fit_decline`, `fit_decline_cum`, `fit_ratio`, `eur`, `forecast` do not. Mixed convention — either all four entry families accept flexible input or none do. Related to API-UX issue #22.

### L8. `ForecastResult.secondary` typed as `Optional[dict]`
No nested typing for `dict[str, dict[str, np.ndarray]]`. Minor — could add a `SecondaryPhase` dataclass per API-UX issue #16.

### L9. `np.array([])` in `field(default_factory=lambda: ...)` is fine but wasted
Lines 96, 153 — empty array default for residuals is allocated on every instantiation. `None` default would be cheaper and semantically clearer ("no residuals yet"). Truly minor.

---

## Summary

21 findings: **2 High, 5 Medium, 9 Low**. Top 3:

1. **H1 — Fitter duplication**: six `_fit_*` helpers share identical R²/residuals boilerplate; two hyperbolic grid-searches share identical scaffold; `fit_decline` + `fit_decline_cum` share identical windowing/dispatch. ~90 lines of duplication total, easily refactored.
2. **M1 — Zero named constants**: 17 magic numbers (b-grid, Duong bounds, curve_fit maxfev, clamp floors) inline in function bodies, directly violating the CLAUDE.md "Named Constants" rule for new code.
3. **M3 + M4 — Sentinel 0.0 & no type hints**: `DeclineResult` / `RatioResult` use 0.0 for method-inapplicable fields (conflates with physical zero); zero type annotations on any function signature across the module.
