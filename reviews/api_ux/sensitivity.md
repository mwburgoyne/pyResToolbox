# API/UX — sensitivity

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/sensitivity/sensitivity.py`, `/home/mark/projects/pyResToolbox/pyrestoolbox/sensitivity/__init__.py`, `/home/mark/projects/pyResToolbox/pyrestoolbox/docs/sensitivity.rst`.

Module is small and clean: 2 functions, 3 dataclasses, 203 lines. `__all__`, re-exports, numpydoc docstrings, and RST are all aligned. The only meaningful gaps are error-message clarity and a handful of missing defensive checks that would make failure modes readable.

## Critical
None.

## High

### 1. `ranges` tuple order is undocumented and unvalidated
`sensitivity.py:151-203`, `sensitivity.rst:101-103`. The RST says "Mapping of parameter names to `(low, high)` tuples" but nothing in the code enforces ordering. A caller passing `{'p': (4000, 1000)}` gets `low_value=4000, high_value=1000` and the sensitivity magnitude is still correct (uses `abs`), but the stored results are semantically backwards. Either document that order is not enforced (and that `low_value` just means "first tuple element"), or validate `lo <= hi`. Related correctness-review item exists; this is the API/UX side.

### 2. `result_key` error handling produces cryptic exceptions
`sensitivity.py:106-112`. If the user passes `result_key='Bo'` and the function returns a dict without that key, the raised `KeyError` comes from the generic dict lookup with no context about which function call failed or at which parameter value. Should be wrapped with a clear "`result_key='Bo'` not found in result for vary_param='p' at value=2000" message.

## Medium

### 1. `sweep`: only `vary_param` presence is checked; `vary_values` is not
`sensitivity.py:138-146`. Guard: `if vary_param not in base_kwargs`. No check that `vary_values` is non-empty, is iterable, or contains values of a reasonable type. An empty list silently yields an empty `SweepResult`. A `None` value yields `TypeError: 'NoneType' is not iterable` with no helpful context.

### 2. `tornado`: empty `ranges={}` silently returns empty entries
`sensitivity.py:176-201`. `base_result` is still computed (one call), and `entries` is an empty list. The caller gets a valid but useless `TornadoResult`. Worth a warning, or an explicit `ValueError("`ranges` must contain at least one parameter")`.

### 3. Type hints mix `list` (lowercase) with `List[TornadoEntry]`
`sensitivity.py:54, 60, 103`. `SweepResult.values: list` and `results: list` use the builtin; `TornadoResult.entries: List[TornadoEntry]` uses `typing.List`. Inside a single module, pick one. For Python 3.9+ (and the package supports cp38+) both work, but inconsistency hurts readability.

### 4. Function signatures type `base_kwargs` as `Dict[str, Any]` but `ranges` as `Dict[str, tuple]` — untyped tuple is weak
`sensitivity.py:115, 151`. `Dict[str, Tuple[float, float]]` would be more honest — callers do pass numeric lo/hi values and there is no generic-object path.

### 5. `sweep` does not return the base result; `tornado` does
Asymmetry. For some sensitivity workflows (e.g., plotting normalized deviation from base) `sweep` users have to separately evaluate `func(**base_kwargs)`. Not broken, but an opportunity for a more consistent API.

### 6. No exception handling in sweep loop
`sensitivity.py:142-146`. If one parameter value raises (e.g., gas Z-factor at `p=-1`), the whole sweep aborts with no partial results. This is a pure UX choice — catching/filling NaN would let users plot partial curves. Covered in correctness review too; repeating because the UX impact is visible.

### 7. Docstring does not explain what happens when `base_result == 0`
`sensitivity.py:189-192`. When `|base_result| == 0` the sensitivity metric silently switches from relative to absolute. This is load-bearing behavior and should be in the docstring/RST so users don't misinterpret tornado bars.

### 8. `SweepResult.values` docstring says "Values that were assigned to the parameter" — RST says "Parameter values used"
`sensitivity.py:56` vs `sensitivity.rst:51-52`. Both fine, but the dataclass stores `list(vary_values)` — if user passed a numpy array, it is coerced to list silently. Document the coercion.

## Low

### 1. `TornadoEntry` has no `__post_init__` validation
`sensitivity.py:65-88`. Dataclass accepts anything — e.g., `low_value=5, high_value=2, sensitivity=-0.3` (negative sensitivity). Not a problem via the public API (the function constructs them correctly), but if users build entries by hand for mocked tests there is no guardrail.

### 2. `sensitivity` absolute-relative formula is only explained in the tornado docstring
`sensitivity.py:170`, `sensitivity.rst:86`. Not in `TornadoEntry.sensitivity` attribute docstring, which says "Absolute relative sensitivity |high_result - low_result| / |base_result|" — that is correct but only visible if you introspect the dataclass. Cross-reference helps.

### 3. `tornado` casts all results to `float` via `float(...)`
`sensitivity.py:174, 183, 187`. If `result_key` returns a non-scalar (list of SegmentResult, array), the cast fails with a naked `TypeError`. Error message should name the parameter and suggest setting `result_key`.

### 4. `__init__.py` / `__all__` — clean
`sensitivity/__init__.py` is a single-line `from .sensitivity import *`. `sensitivity.py:37-40` declares `__all__` with all 5 public names. Good.

### 5. Neither `sweep` nor `tornado` documents thread-safety
Concern is low (pure function evaluation), but if `func` has internal state users should know the loop is sequential. One-liner in the docstring.

### 6. RST examples depend on `gas.gas_z` values matching exactly
`sensitivity.rst:71-74, 130-140`. Hardcoded `0.9260188251531628` and `0.0886413428394408`. These are captured via `test_doc_examples.py`. Any change to `gas_z` defaults invalidates the doc — this is by design per CLAUDE.md, just worth flagging that sensitivity's doc examples are coupled to gas frozen baselines.
