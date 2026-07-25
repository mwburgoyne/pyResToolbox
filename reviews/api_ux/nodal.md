# API/UX — nodal

Scope: review brief listed `fbhp, fthp, vlp, gradient, operating_point` — but only `fbhp, outflow_curve, ipr_curve, operating_point` exist (plus `WellSegment, Completion, Reservoir, NodalResult`). Absence of `fthp` (reverse solve for wellhead pressure) and any stand-alone `gradient`/`vlp` helper is itself flagged below.

## Critical

### 1. Invalid-method error message doesn't list valid options
`validate_methods()` (`pyrestoolbox/validate/validate.py:28`) raises:
`ValueError("An incorrect method was specified: 'INVALID' for 'vlpmethod'")`
No hint that the accepted tokens are HB / WG / GRAY / BB. For an LLM agent or first-time user this is the single worst UX in the module — the only way to discover valid methods is to read the RST. Include the valid set in the message: `"Invalid vlpmethod 'INVALID'. Valid options: HB, WG, GRAY, BB"`.

### 2. `well_type` is never validated
`fbhp`, `outflow_curve`, `operating_point` accept `well_type='gas'` as default and branch on `== 'gas'` / `== 'oil'`. `ipr_curve` additionally accepts `'water'`. A typo (`'Oil'`, `'OIL'`, `'liquid'`, `'condensate'`) silently takes the `oil`/default branch and returns plausible-looking numbers. No `ValueError`. This is a silent wrong-answer path. Validate once at entry with an explicit whitelist.

## High

### 3. No `fthp` (reverse direction) function
There is no public "solve for tubing head pressure given BHP" (the natural inverse of `fbhp`). Field work routinely needs this — wellhead choke sizing, back-calculation of surface pressure from DHSG data, matching against shut-in WHP. Users are left to wrap `bisect_solve` themselves around `fbhp`. Either add `fthp(bhp, completion, ...)` or document clearly that only forward (THP → BHP) solves are supported and explain the `injection=True` workaround for downward marches.

### 4. No pressure-traverse / gradient output
`fbhp` returns only the bottom-hole scalar. There is no way to recover the depth-vs-pressure traverse, hydrostatic/friction split, or per-segment holdup without forking the internal segment march. This is the most-requested VLP diagnostic output in practice and is trivial to expose (the segment loop already computes everything). Suggest a `return_profile=False` kwarg that returns `NodalResult` with `md`, `tvd`, `p`, `dpdz_hydro`, `dpdz_fric`, `holdup` arrays. Related: no `gradient(completion, ...)` helper despite being listed in brief.

### 5. No survey-format wellbore input
A deviated completion must be constructed by manually assembling `[WellSegment(md=..., deviation=...), ...]` tuples where each `md` is *segment length*, not *cumulative MD*. Every field survey in existence is `(MD, inclination)` pairs with MD cumulative. RST does not document this convention clearly — "Measured depth of this segment" could plausibly read as cumulative. Add a `Completion.from_survey([(md1,inc1),(md2,inc2)...], tid, tht, bht)` classmethod that converts cumulative MD → segment lengths and piecewise-constant inclination.

### 6. `rates` parameter in `outflow_curve` is ambiguous under `metric=True`
Docstring: `"rates: Explicit list of rates to evaluate (same units as output rates)"`. Code (line 1837): `rate_mmscfd = rate * SM3_TO_MMSCF if metric else rate`. So a metric user passing `rates=[50000, 100000]` intending sm3/d works, but there is zero validation or documentation that these are sm3/d in metric mode — easy to confuse with MMscf/d. Adding a doc-example in metric mode would help, or reject magnitude outliers.

### 7. Gas IPR in Mscf/d vs VLP in MMscf/d is a persistent footgun
RST has a warning box calling this out, `operating_point` handles the 1000× scaling internally — so the library is aware. But `ipr_curve()` and `outflow_curve()` still return incompatible units by default. If a user plots both dicts directly, their VLP curve looks flat at rate < 0.05 against a 13,000 Mscf/d IPR axis. Recommend either (a) unify to MMscf/d in both, or (b) add an explicit `rate_units='MMscf/d'` kwarg to `ipr_curve` for gas, or (c) return unit labels in the dict.

### 8. `api` has wrong default for oil wells
Default `api=45` is appropriate for condensate (gas well default) but unusually light for an oil well. On `fbhp(..., well_type='oil', ...)` without `api` specified, users silently get a 45-API oil — easy to miss since no warning fires. Either (a) make `api` required (no default) when `well_type='oil'` and no `oil_pvt` supplied, or (b) raise a warning when API > 40 and well_type='oil'.

### 9. Error messages omit unit context
Example: `WellSegment(md=-1, id=2.441, metric=True)` raises `"Measured depth md must be positive, got -0.3048"` — the user sees the post-conversion ft value, not the -1 m they typed. All WellSegment/Completion/Reservoir validators echo the *stored* (oilfield) value regardless of which unit system the caller used. Fix: validate before conversion, or echo both (`got -1 m`).

### 10. `WellSegment(id=...)` shadows Python builtin
`id` is a Python builtin; using it as a positional arg name means IDEs / linters may flag every call. Low user-impact but a common style nit. Consider `diameter=` or `internal_diameter=` as the canonical name with `id=` kept as a deprecated alias.

## Medium

### 11. `oil_pvt` is accepted but `gas_pvt` is unused by VLP
`fbhp` docstring: `"gas_pvt: GasPVT object (unused by VLP methods directly, reserved for future use)"`. Accepting a parameter that does nothing is worse than not accepting it — users pass `gas_pvt` thinking it propagates composition (CO2/H2S/N2/H2) into the VLP calculation and silently get a `gsg`-only result. Either wire it through (the internal `_z_factor` helper uses `_sutton_tc_pc(sg)` with no inert correction) or drop the parameter.

### 12. Duplicated & overlapping parameters when `oil_pvt` is passed
`fbhp(oil_pvt=opvt, api=35, sgsp=0.65, pb=2500, rsb=500)` — caller passes both `oil_pvt` (which contains api/sgsp/pb/rsb) AND scalar overrides. Code silently overwrites the scalars from `oil_pvt.*` (line 1714-1719). No warning on conflict. User intent is ambiguous. Either (a) warn if both are specified and differ, or (b) document that `oil_pvt` wins.

### 13. `operating_point` silently returns `op_rate=0` on bisection failure
Line 2070-2073:
```python
try:
    op_rate = bisect_solve(...)
except (RuntimeError, ValueError):
    op_rate = 0.0
```
When no operating point exists (VLP everywhere above IPR), the function returns `{'rate': 0.0, 'bhp': pr, ...}` with no indication. User gets the reservoir pressure as BHP and zero rate — looks plausible but tells them nothing. Warn or add a `converged: bool` field in the result.

### 14. `injection=True` is an afterthought, undocumented workflow
It's a boolean kwarg on `fbhp`, flipping a sign inside the gradient. But an injector is a very different workflow (water/gas injection, different PVT, THP > BHP means pressure *rising* downward). No dedicated injection example in RST. `outflow_curve` and `operating_point` don't forward `injection` to the inner `fbhp` calls (lines 1841, 1847, 2048, 2054, 2077, 2082) — those are permanently `injection=False`. So `injection=True` is only usable via direct `fbhp` calls. Either wire it through or document the restriction.

### 15. `oil_vis` and `uo` are duplicate concepts
`fbhp(oil_vis=1.0)` — condensate viscosity, only used for gas wells. `ipr_curve(uo=1.0)` — oil viscosity, only used for oil wells. `operating_point` has *both*. Same quantity, different names, different well types. Confusing. Standardize on one name with well-type context.

### 16. `Completion` has two construction modes with surprising attribute semantics in segment mode
Line 385-390: in segment mode, `self.tid = self._segments[0].id` and `self.length = self._segments[0].md`. A 3-segment multi-diameter completion sets `tid` and `length` from segment 0 only — users may read `completion.tid` and be misled. These "legacy attributes for compatibility" should either be removed or raise on access in segment mode.

### 17. `max_rate` default of 50 MMscf/d and 10000 STB/d is hard-coded
`outflow_curve` auto-generates rates 0.01–50 MMscf/d (gas) or 1–10000 STB/d (oil). No warning if the operating point is outside this range. Small wells (< 1 MMscf/d) get only one useful rate point in the 20-point default. Recommend: compute a sensible max from AOF of associated IPR, or warn when all BHP values exceed reservoir pressure.

### 18. `outflow_curve` returns `'rates'` (plural), `ipr_curve` returns `'rate'` (singular)
Minor inconsistency in dict keys. `operating_point` returns `'rate'` (scalar, singular). A user post-processing both dicts has to remember which is which. Standardize.

### 19. `pr=0 disables` condensate dropout — magic sentinel
Docstring: `"pr: Reservoir pressure for condensate dropout. 0 disables"`. Using 0 as a magic-value for "off" is fragile (a user with a depleted tight-gas well at 50 psi would need to handwave). Use `None` for "disabled" and make 0 a valid-but-meaningless value, or add a `condensate_dropout=True/False` kwarg.

## Low

### 20. `NodalResult` supports attribute access but `vlp`/`ipr` nested dicts don't
`result.rate` works (line 83-87). But `result.vlp` returns a plain dict, so `result.vlp.rates` raises. Inconsistent dot-access across nesting levels. Wrap nested dicts in `NodalResult` too.

### 21. `Reservoir` accepts `D` (non-Darcy) only in oilfield units per docstring, but metric conversion exists
Line 582-583: `D` is converted with `D_PER_SM3_TO_D_PER_MSCF` when `metric=True`. Docstring *does* mention metric D units. Fine — but the RST table says `"day/mscf for gas, or day/sm3 if metric=True"` with no mention that `D` for oil (where it should be zero) still silently converts. Low-impact clarification.

### 22. No docstring example for `outflow_curve` or `ipr_curve` metric mode
RST gives field-unit examples only. Metric users must infer. Add one metric example for each.

### 23. `n_rates` vs `n_points`: different argument names for the same concept
`outflow_curve(n_rates=20)`. `ipr_curve(n_points=20)`. `operating_point(n_points=25)`. Three functions, two names, same concept (size of the sampled curve). Unify.

### 24. Unit tables don't specify pressure convention at tubing shoe vs wellhead
All examples use `thp` (tubing head) but the term "bottom hole" in `fbhp` is ambiguous between tubing shoe and mid-perforation when `mpd > length`. Completion class docstring covers this but `fbhp` docstring does not reference `mpd`. Add one-line note.

### 25. Default `vlpmethod='WG'` — OK but unmentioned rationale
RST default-picks WG (drift-flux, works for all deviations) over HB (vertical-only). Fine choice but no justification in the docstring. Add one-line: `"default 'WG' because it handles all inclinations"`.

### 26. `pr` and `reservoir.pr` can disagree silently
`fbhp(..., pr=3500, ...)` (for condensate dropout) can be called independently of the `Reservoir` object's `pr`. `operating_point` forwards `reservoir.pr` into `fbhp(pr=pr)` (line 2011, 2050) — good. But a user calling `fbhp` directly with `pr=3500` and also passing a `Reservoir(pr=3000, ...)` somewhere else in their workflow has two separate pr values. Accept only one. (Caveat: `fbhp` doesn't take `reservoir`, so cross-contamination is narrow — still worth a doc note.)

### 27. `NodalResult` docstring lists three factory functions but `fbhp` returns a plain float
Line 76-82: "Dict subclass returned by outflow_curve, ipr_curve, and operating_point." `fbhp` returns float. Fine, but this pattern difference (scalar vs dict return) should be called out as "all nodal functions return NodalResult *except* `fbhp` which returns a single scalar BHP."

---

## Summary
27 issues: 2 Critical, 8 High, 9 Medium, 8 Low. Top 3:
1. Invalid-method error doesn't enumerate valid options — worst agent-facing UX in the module.
2. `well_type` never validated — typos silently dispatch to default branch.
3. No pressure-traverse / gradient output and no `fthp` reverse solve — two most-requested VLP workflows absent despite internals being trivial to expose.
