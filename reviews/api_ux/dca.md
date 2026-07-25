# API/UX — dca

## Critical

### 1. Inconsistent casing in arg names: `qi`, `di`, `Np`, `Qcum`, `t_calendar`, `b`, `q_min`
`fit_decline_cum(Np, ...)` uses PascalCase for the array argument while `fit_decline(t, q, ...)` uses lowercase. `ForecastResult.Qcum` also breaks PEP 8 / convention used elsewhere in the package (`gas.gas_z`, etc. all lowercase). Jarring for API consumers and autocomplete. Either standardize on lowercase (`np_`, `qcum` — though `np` collides with numpy) or document explicitly why engineering-notation casing is preferred. Current mix is the worst of both.

### 2. `forecast()` silently starts the time grid at `t=dt`, not `t=0`
`t = np.arange(dt, t_end + dt/2, dt)` — no data point is produced at t=0, and the user has no way to override this. This means `Qcum[0] = q(dt)*dt`, not zero, and `fc.t[0]` is never zero. Not documented. Users comparing `fc.Qcum` to `arps_cum(qi, di, b, fc.t)` get off-by-one-step results. Either document prominently or start at t=0.

## High

### 3. Time-unit ambiguity is noted in one place but invisible at call-site
`dca.rst` says "unit-agnostic … if you pass stb/d and months, qi comes back in stb/d and di in 1/month". Good. But individual docstrings for `arps_rate`, `duong_rate`, `eur`, `forecast` only say "(1/time)" or "Time". A user reading only a docstring (e.g. via LSP hover) can't tell whether `di=0.1` means 0.1/day or 0.1/month. At minimum, every public docstring should include a one-line "unit-agnostic" reminder.

### 4. `forecast()` does not accept a `DataFrame` output option, and no `arrays-in` variant
Input to `fit_decline` / `fit_decline_cum` is strictly `(t, q)` arrays — no DataFrame convenience (e.g., `fit_decline(df, t_col='date', q_col='oil_rate')`). `ForecastResult` also returns raw arrays, not a DataFrame. Petroleum DCA workflows are overwhelmingly DataFrame-based (monthly production). Forcing users to unpack columns and repack results is friction. Either add `as_dataframe()` helpers on result classes, or document that this is intentionally array-only.

### 5. `uptime` parameter has ambiguous semantics in `forecast()` vs `fit_decline_cum()`
`forecast(uptime=0.9)` = scalar applied uniformly. `fit_decline_cum` returns `uptime_mean` and `uptime_history` (per-interval). There's no way to feed `uptime_history` back into `forecast()` as a time-varying schedule. User has two parallel concepts of uptime that don't connect. Either accept an array/callable for `forecast(uptime=...)` or document the limitation loudly.

### 6. Result dataclasses have sentinel `0.0` for unused fields
`DeclineResult.a=0.0, m=0.0` when method is Arps; `.di=0.0, .b=0.0` when method is Duong; `RatioResult.c=0.0` unless logistic. Looks like a real parameter value but is a sentinel. Prefer `None` for unused fields, or use method-specific subclasses. Currently users must know the method to interpret the fields — easy to misread `di=0.0` as "zero decline" for a Duong result.

### 7. `RatioResult.a` / `.b` / `.c` are positional overloaded with method-dependent meanings
Docstring: `a = Primary parameter (intercept / coefficient / Rmax for logistic)`. Same field means four different things depending on `method`. This is a footgun — anyone post-processing fits will write `if rr.method == 'linear': gor_intercept = rr.a`. Much safer to expose named properties (`.intercept`, `.slope`, `.Rmax`, `.rate_constant`) or a `params: dict` attribute.

## Medium

### 8. `fit_decline(method='best')` gives no tie-break or reporting
Returns single best-R² result; the user cannot see which other methods were tried, their R², or how close they were. For `best` fits, a `candidates: list[DeclineResult]` attribute on the winner would aid trust and diagnostics. Also — no warning when the "best" R² is poor (e.g., < 0.5).

### 9. `duong_rate` argument order `(qi, a, m, t)` is inconsistent with literature
Duong papers conventionally write `q = q1 * t^(-m) * exp(...)` with `m` as the first shape parameter. Order `(a, m)` here is fine but `a` represents the intercept of a log-log plot with different conventions across SPE papers. At a minimum, clarify in docstring what definition of `a` and `m` is used (e.g., "Duong 2011, SPE 137748, Eq. X").

### 10. `eur()` returns scalar — but `q_min` units must match `qi`, unchecked
Signature gives no hint that units must match. `dca.eur(qi=1000, di=0.1, b=0, q_min=10)` — is 10 an economic limit in bbl/d, bbl/month, or Mscf/d? The function has no way to know. Compounded by time-unit ambiguity (#3). Docstring should warn.

### 11. `fit_ratio` method names inconsistent with `fit_decline`
`fit_decline`: `{exponential, harmonic, hyperbolic, duong, best}`. `fit_ratio`: `{linear, exponential, power, logistic, best}`. Both have `exponential` but with different meanings (decline vs growth). Consider namespaced methods or document prominently.

### 12. `ForecastResult.eur` is "cumulative at end of forecast", not true EUR
The standalone `eur()` function returns EUR at economic limit. `ForecastResult.eur` is just `Qcum[-1]`, which may or may not represent EUR depending on whether `q_min` triggered truncation. Confusing field name. Suggest `cumulative_final` or document the distinction in the docstring (RST does say "final cumulative value" — good — but the dataclass docstring says "Estimated ultimate recovery" which is the misleading one).

### 13. Error messages: some are actionable, others not
Good: `"b must be between 0 and 1, got {b}"`. Bad: `"All decline fitting methods failed"` — no hint at why (too few points? all zero? monotonic?). Bad: `"Cumulative decline fitting with method '{method}' failed"` — should say what to try next (different window? different method?).

### 14. `fit_decline` window shift is documented but surprising
`t_start=20, t_end=60` silently shifts `t -> t - t[0]`. Docstring in code says "# Shift so window starts at t=0"; RST says "returned qi represents the rate at the window start". This means `result.qi` is NOT the rate at `t=0` of your original data — it's the rate at `t=t_start`. Easy to get wrong. Suggest explicit parameter `shift_to_origin=True` so behavior is opt-in and obvious.

## Low

### 15. `ratios` dict in `forecast()` is stringly-typed
`ratios={'GOR': rr}` — the key becomes a column name in `secondary`. No validation that the names are sensible. A small enum or typed container would help, but this is minor.

### 16. `ForecastResult.secondary` nested dict shape undocumented in type hints
`Optional[dict]` is the only type hint. The actual shape `dict[str, dict[str, np.ndarray]]` with keys `'ratio'/'rate'/'cum'` is documented in the docstring but not typed. Consider a `SecondaryPhase` dataclass.

### 17. `arps_cum(qi=..., di=..., b=..., t=...)` is keyword-friendly, but Arps b-factor is often called `b` in decline papers and `n` in others
Minor stylistic issue; current naming is fine but a b-vs-n comment in docstring would help newcomers.

### 18. No `__repr__` overrides on dataclasses — arrays print fully
`print(result)` on a `DeclineResult` with a 500-element residual array floods the terminal. Default dataclass repr includes full array. Add `repr=False` on `residuals` / `uptime_history` fields, or custom `__repr__`.

### 19. `fit_decline_cum` uptime inference silently returns `uptime=1.0` on `dt_cal <= 0`
If calendar times are duplicated or decreasing, each such interval gets `uptime=1.0` (perfect) rather than a warning or skip. Silent cleansing of bad input. Minor but could surprise.

### 20. `ratio_forecast(rr, 10.0)` returns scalar 6.0 (doc example) — fine, but no clipping on logistic/power domain
If you call `ratio_forecast` with `x < 0` on a `power` model, you get `nan` or complex — no friendly error. Edge case but worth a guard.

### 21. `forecast(result, t_end, dt=1.0, ...)` — `t_end` unit is same as `result.di`'s implied time unit, but there's no cross-check
If you fit with `di=0.05` (per month) and call `forecast(result, t_end=50, dt=1.0)` thinking `t_end=50 years`, you'll get 50 months. No way for the library to catch this. Docstring mention only.

### 22. `arps_rate` uses `convert_to_numpy/process_output` but `fit_decline` does not accept scalars
Asymmetric. `arps_rate` accepts list/array/scalar; `fit_decline` requires array-like length ≥ 3. Fine functionally, but inconsistent API shape across sibling functions.

---

## Summary
21 issues: 2 Critical, 5 High, 7 Medium, 7 Low. Top 3:
1. Inconsistent arg casing (`Np`, `Qcum` vs `t`, `q`, `qi`) — pervasive and user-facing.
2. `forecast()` skips t=0, producing subtle off-by-one-step mismatches with analytical `arps_cum()`.
3. Time-unit ambiguity — unit-agnostic design is fine, but individual docstrings don't remind users, and `eur()` / `forecast()` inputs have no sanity check.
