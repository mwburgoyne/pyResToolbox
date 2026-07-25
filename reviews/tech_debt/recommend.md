# Tech Debt — recommend

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/recommend/recommend.py` (275 lines), `/home/mark/projects/pyResToolbox/pyrestoolbox/tests/test_recommend.py` (13 tests).

Small, self-contained module: 4 functions + 1 dataclass, no external deps beyond `dataclasses` and `typing`. Main tech-debt items are unused parameters (carried for API-compatibility but never consumed), hand-curated alternatives lists that duplicate info already in `pyrestoolbox.classes` Enums, and magic thresholds.

## High

### 1. `sg` parameter is dead — never used in any decision branch
`recommend.py:71, 237`. Declared in both `recommend_gas_methods` and `recommend_methods`, documented, defaulted to `0.65`. No branch in `recommend_gas_methods` body (lines 93-149) reads `sg`. Threaded through `recommend_methods` unchanged. Same pattern flagged in API/UX review; from tech-debt angle it is just unused-parameter noise that has to be kept in the docstring + signature in sync forever.

Either wire `sg` into the Z-factor recommendation (e.g., DAK loses accuracy for `sg > 0.85`) or remove the parameter and delete the threading through `recommend_methods`.

### 2. `well_type` parameter is dead — same problem
`recommend.py:199, 241, 273`. `recommend_vlp_method` accepts `well_type` and `recommend_methods` threads it through, but the function body at lines 214-232 does not reference it. Same correction options as #1.

### 3. `alternatives` string lists are hand-curated and drift-prone against `classes.py` Enums
`recommend.py:99, 116, 121, 129, 134, 142, 147, 181, 186, 192, 223, 231`. Hard-coded strings like `'DAK'`, `'PMC'`, `'BNS'`, `'VELAR'`, `'WG'`, `'HB'`, `'GRAY'`, `'BB'` — these must match Enum member names in `pyrestoolbox.classes`. Adding a new method to `classes.z_method` requires a manual sweep of `recommend.py` to decide whether it becomes an alternative in any branch.

**Fix:** import Enum members and reference `z_method.DAK.name` etc. This at least breaks loudly (AttributeError) instead of silently omitting the method if an Enum is renamed.

## Medium

### 1. Magic thresholds `0.55`, `0.10`, `30` (degrees) are inlined
`recommend.py:111, 124, 216`. CLAUDE.md's "Named Constants" rule applies. Extract:
```python
# Inert fraction above which BNS 5-component PR-EOS is preferred over DAK/PMC
_INERT_THRESHOLD_BNS = 0.55
# CO2 or H2S fraction above which PMC critical-property corrections are material
_IMPURITY_THRESHOLD_PMC = 0.10
# Wellbore deviation (deg) above which HB and Gray become unreliable
_VLP_DEVIATION_LIMIT = 30.0
```
Cite the source (Burgoyne/BNS paper, industry rule of thumb, SPE monograph) inline.

### 2. Rationale-string duplication — same compositional branch logic spelled twice
`recommend.py:97-110` (H2 branch zmethod+cmethod) and `111-123` (high-inert branch zmethod+cmethod) each build two separate `MethodRecommendation` objects with mostly-parallel rationale strings. Five branches × 2 dicts = 10 near-duplicated constructor calls. A small factory like
```python
def _make_pair(z_rec, z_rationale, c_rec, c_rationale, ...):
    return {'zmethod': ..., 'cmethod': ...}
```
would halve the repetition and make branch intent cleaner.

### 3. Test coverage gaps
`test_recommend.py` has 13 tests. Covered: clean gas, H2, high inerts, moderate CO2, moderate H2S, oil default/heavy/light, VLP vertical/deviated/30deg, master gas-only, master with-oil. Gaps:
- `recommend_gas_methods` with N2 only (no CO2/H2S) — untested
- `recommend_gas_methods` with `sg=0.9` (heavy gas) — untested (but parameter is dead anyway, see High #1)
- `recommend_oil_methods` with mid-range API (e.g., 30) vs heavy vs light — tests only check the `recommended` string, not the `rationale` differentiating text
- `recommend_vlp_method` with `well_type='oil'` — untested (dead param, see High #2)
- `recommend_methods` with no args at all — default-path smoke test missing
- `MethodRecommendation` dataclass direct construction / equality / `mandatory` flag — untested
- No error-case tests (negative fractions, etc.) — but the function has no validation to test, so this is latent tech-debt

### 4. `mandatory` flag semantics inconsistently documented
`recommend.py:48-68`. The dataclass docstring says `mandatory: bool — True if the recommended method is the only valid choice.` But in the H2 branch (line 102) `alternatives=[]` is paired with `mandatory=True` — so the flag is technically redundant with "empty alternatives". If `mandatory=True` is the contract, testing `len(alternatives) == 0` should be equivalent. Pick one source of truth.

### 5. Return-dict shape varies by input
`recommend.py:267-275`. `recommend_methods` returns 3 keys without `api`, 6 keys with. Dict-shape polymorphism. Either always return 6 keys (missing ones as `None`) or split into two functions. Callers doing `recs['pbmethod']` without checking raise `KeyError` — tech-debt in the sense that downstream code has to defensively check.

## Low

### 1. `sensitivity.py` and `recommend.py` both use `typing.List` / `typing.Dict` / `typing.Optional`
`recommend.py:44`. Fine for cp38, but Python 3.9+ allows `list[str]`, `dict[str, MethodRecommendation]`, `str | None`. Codebase targets cp38-cp313 (per CLAUDE.md) so stick with `typing.*` — but be aware if support floor moves to 3.10.

### 2. `field(default_factory=list)` at line 67 is correct pattern — good
Noted for completeness.

### 3. Percentage formatting `{h2:.1%}` in rationale strings will show `0.0%` for zero values — but only reachable if the branch ever allows it
`recommend.py:100, 115, 128`. H2 branch gates on `h2 > 0` (good). Inerts branch doesn't gate individual components but the rationale at line 115 formats only `inerts`, not individual `co2/h2s/n2/h2`. Currently safe; noted as a potential footgun if rationale strings are expanded.

### 4. Module-level `__all__` at lines 37-41 matches the documented Functions/Classes section at lines 26-34 — good
No drift.

### 5. All four functions return `Dict[str, MethodRecommendation]` but the dict keys are unconstrained (bare strings)
Could be strengthened with a `Literal['zmethod', 'cmethod', ...]` or `TypedDict` in a future refactor. Low priority — keeps the API flexible today.

### 6. `MethodRecommendation.alternatives` order is undocumented
Tech-debt adjacent to API/UX — if the order is priority-ranked (the code appears to intend that), add one-liner in the dataclass docstring. Prevents future contributors adding alphabetical-sort calls and breaking implicit contract.

### 7. No module-level import of `classes` Enums
`recommend.py` touches Z-factor/oil/VLP methods without importing from `pyrestoolbox.classes`, which means it cannot validate that `'DAK'` is actually a valid `z_method` name at definition time. See High #3.
