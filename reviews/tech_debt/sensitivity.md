# Tech Debt — sensitivity

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/sensitivity/sensitivity.py` (202 lines), `/home/mark/projects/pyResToolbox/pyrestoolbox/tests/test_sensitivity.py` (12 tests).

**This module is clean.** 2 functions + 3 dataclasses + 1 private helper, no external dependencies beyond `dataclasses` and `typing`, minimal duplication, no magic-number problems. All reasonable tech-debt issues are minor or covered elsewhere. Nothing High.

## High

None.

## Medium

### 1. Duplicated per-parameter call-site construction in `tornado()`
`sensitivity.py:181-187`:
```python
kwargs_lo = dict(base_kwargs)
kwargs_lo[param] = lo
lo_result = float(_extract_result(func(**kwargs_lo), result_key))

kwargs_hi = dict(base_kwargs)
kwargs_hi[param] = hi
hi_result = float(_extract_result(func(**kwargs_hi), result_key))
```
Same three-line pattern for lo and hi. Extract a local helper:
```python
def _eval_at(param, value):
    kwargs = dict(base_kwargs)
    kwargs[param] = value
    return float(_extract_result(func(**kwargs), result_key))
```
Then both call sites become `lo_result = _eval_at(param, lo)`. Six lines collapse to two; the intent ("evaluate func with one override") becomes visible. Same pattern exists in `sweep` (lines 143-146) — the helper could be reused.

### 2. Mixed `list` vs `List[...]` type hints inside the same file
`sensitivity.py:54, 60, 103`. `SweepResult.values: list` and `SweepResult.results: list` use the built-in bare `list`. `TornadoResult.entries: List[TornadoEntry]` uses `typing.List`. Pick one. (Codebase supports cp38+, so `typing.List` remains compatible — use that for consistency with `recommend.py`, which also uses `typing.List`.)

### 3. `ranges: Dict[str, tuple]` — untyped tuple
`sensitivity.py:152`. Should be `Dict[str, Tuple[float, float]]` for agent/IDE type-assistance. Caller always passes `(float, float)`; the bare `tuple` hint loses that.

### 4. Test coverage gaps
`test_sensitivity.py` has 12 tests. Gaps:
- Neither `sweep` nor `tornado` tested with a function that returns an **object with attribute access** (`_extract_result` line 112 uses `getattr`) — only dict-return tested.
- `sweep` with `result_key` pointing at missing key — error path uncovered.
- `tornado` with `ranges={}` (empty dict) — result shape uncovered.
- `tornado` sorting with three or more entries where the middle one has smallest sensitivity — only 2-entry ordering is tested.
- `SweepResult`/`TornadoEntry`/`TornadoResult` dataclass equality, repr, `field(default_factory=list)` default — untested.
- `_extract_result` private helper not directly tested — fine, but the `if result_key is None` branch is indirectly covered while the `getattr` branch is not.

### 5. `_extract_result` signature asymmetry vs callers
`sensitivity.py:106-112`. Accepts `(result, result_key)`. Three call sites (lines 146, 173, 183, 187) pass `result_key` through. Internal only — fine — but the helper does not validate that `result_key` is a string when the `result` is a dict/object; a None result silently raises `TypeError`. Either document or validate.

## Low

### 1. `result_key: Optional[str] = None` — but `_extract_result` doesn't enforce string
`sensitivity.py:117, 153`. A list of keys (`result_key=['a', 'b']`) would pass the `None` check and then fail deep inside `result[result_key]` with a cryptic TypeError. Since neither multi-key extraction nor nested-path extraction is supported, this is purely a type-hint-honoured-at-runtime gap. Trivial to tighten.

### 2. `sweep()` does not return the base result; `tornado()` does
`sensitivity.py:148, 202`. Asymmetric API between the two. Not a tech-debt blocker but a maintenance itch — if a future caller wants "sweep result deltas from base", they have to compute base separately. Flagged also in API/UX review.

### 3. `list(vary_values)` at line 148 silently coerces numpy arrays
`sensitivity.py:148`. `np.ndarray` → Python list is fine but loses dtype and may surprise users plotting `result.values` as a numpy array later. Document or use `np.asarray` and store raw.

### 4. `TornadoEntry` has 6 fields; dataclass generates default `__repr__` that prints on one line
`sensitivity.py:65-88`. For interactive use (`print(tornado(...))`) the single-line repr of `TornadoResult` containing N entries is hard to read. Not tech-debt in the strict sense — cosmetic.

### 5. No shared helper between `sweep` and `tornado`
Both build `kwargs = dict(base_kwargs); kwargs[name] = val; _extract_result(func(**kwargs), result_key)` around a loop. One file, one helper — worth factoring out once Medium #1 lands.

### 6. `base_kwargs` type hint is `Dict[str, Any]`
`sensitivity.py:115, 151`. Reasonable — user functions accept arbitrary kwargs — but any downstream check against known parameter names would benefit from narrower typing (`Mapping[str, object]` to signal read-only intent; the function shallow-copies on every iteration anyway). Micro-nit.

### 7. `sensitivity/__init__.py` is a single `from .sensitivity import *` — matches the pattern of the other small modules
Noted only to confirm the package wiring is fine. No action.
