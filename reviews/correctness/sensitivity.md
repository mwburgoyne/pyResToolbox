# Correctness — sensitivity module

Reviewed `/home/mark/projects/pyResToolbox/pyrestoolbox/sensitivity/` (203 lines, 2 classes, 2 functions).

## Critical (bug, user gets wrong answer)
_None._

## High (edge case, silent failure)

- **sensitivity.py:189-192, 201** — NaN/Inf `base_result` silently produces unpredictable tornado sort order.
  NaN comparisons fail → `abs(base_result) > 0` is False → falls through to `sensitivity = abs(hi_result - lo_result)` which can itself be NaN. No warning raised.
  **Fix:** check `np.isfinite(base_result)` explicitly; raise or warn if not.

- **sensitivity.py:177-199** — No validation that ranges tuples follow `(low, high)` order.
  If caller passes `(5.0, 2.0)`, `TornadoEntry.low_value=5.0`/`high_value=2.0` — semantically backwards.
  **Fix:** assert `lo <= hi` with a helpful error message, or auto-sort.

## Medium (validity / robustness)

- **sensitivity.py:142-146** — `func(**kwargs)` in sweep loop has no exception handling. One failing point aborts the whole sweep with no partial results. Inconsistent with gas/oil/brine's fault tolerance. **Fix:** catch exceptions, record NaN, emit warning, continue.
- **sensitivity.py:106-112** — `result_key` extraction raises cryptic KeyError/AttributeError if typo'd. **Fix:** wrap with a clear error message naming the missing key.
- **sensitivity.py:174, 183, 187** — `float(...)` coercion has no type check first. TypeError if extraction returns str/list/None.

## Low (nitpick)

- **sensitivity.py:143** — `kwargs = dict(base_kwargs)` is shallow; mutable values shared across sweep points. Rarely hits in practice.
- **sensitivity.py:177** — Empty `ranges={}` silently returns empty TornadoResult. Probably indicates user error; worth a warning.

## Notes — algorithms are correct

Sweep grid generation, sensitivity formula (`|hi-lo|/|base|`), descending sort, cartesian product over multi-param — all correct on the happy path. Issues are all defensive-check gaps.
