# Tech Debt — layer

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/layer/layer.py` (298 lines), `/home/mark/projects/pyResToolbox/pyrestoolbox/layer/__init__.py` (empty), `/home/mark/projects/pyResToolbox/pyrestoolbox/tests/test_layer.py` (12 tests).

Small module, mostly clean internals. The main tech debt is duplication of the EXP/LANG dispatch logic across four public functions and a pile of un-named magic constants (the clamping limits 709, 25000, 0.000001, and the small/large Lorenz short-circuit boundaries).

## High

### 1. EXP/LANG dispatch logic is duplicated four times
The block
```python
B = max(B, 0.000001)
if method == "EXP":
    B = min(B, 709)
    ...  # EXP formula
else:
    B = min(B, 25000)
    PL = 1 / B
    VL = PL + 1
    ...  # LANG formula
```
appears — with small variations — in `lorenz2b` (LorenzErr inner, lines 85-96), `lorenzfromb` (lines 116-124), `lorenz_from_flow_fraction` (BErr inner, lines 165-174), `lorenz_2_flow_frac` (lines 206-214), and `lorenz_2_layers` (lines 258-288). Five copies of the same clamping + two-branch formula.

**Fix:** extract two private helpers — `_lorenz_from_B(B, method)` and `_flow_frac_from_B(B, phih, method)` — that handle clamp + dispatch. Every caller collapses to one call.

### 2. Magic numbers for clamping and near-boundary short-circuits
`layer.py:66-74, 77-80, 87-90, 116-123, 159, 167-170, 206-211, 258-262`. Hard-coded:
- `709` (EXP B upper limit — `np.exp(709)` is near float64 max)
- `25000` (LANG B upper limit)
- `0.000001` (B lower limit)
- `2 / 1000` / `1 / 1000` (returned B for near-homogeneous EXP/LANG)
- `0.000333` (lower Lorenz short-circuit threshold)
- `0.997179125528914` (upper Lorenz short-circuit threshold)

None are named. CLAUDE.md's "Named Constants" rule applies: extract `_B_MAX_EXP = 709`, `_B_MAX_LANG = 25000`, `_B_MIN = 1e-6`, `_LORENZ_HOMOGENEOUS_THRESHOLD = 0.000333`, `_LORENZ_HETEROGENEOUS_THRESHOLD = 0.997179125528914`. Cite the LinkedIn article or a more durable source (see API/UX review).

## Medium

### 1. Unused import `numpy.typing as npt`
`layer.py:40`. Imported but never referenced. Flagged by `pyflakes`/`ruff F401`. Remove.

### 2. Duplicate Langmuir formula in docstrings
The exact block
```
For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
Where PL = 1 / B and VL = PL + 1
```
is copy-pasted verbatim into every function's docstring (lines 51-53, 109-111, 138-140, 192-194, 240-242). Five copies. If the formula or citation ever changes, all five must be edited together. Either hoist to a module-level docstring cross-referenced with a one-line pointer, or leave a single canonical copy in `lorenz2b` and have other docstrings say "see lorenz2b".

### 3. `lorenz_2_layers` uses a manual Python comprehension where numpy would be natural
`layer.py:271-294`:
```python
phih = [0] + [sum(phi_h_fracs[: i + 1]) for i in range(len(phi_h_fracs) - 1)] + [1.0]
```
is `np.concatenate([[0], np.cumsum(phi_h_fracs)[:-1], [1.0]])` — and `np.cumsum` avoids the quadratic repeated `sum(phi_h_fracs[:i+1])`. For the `else` branch at line 278, `np.arange(0, 1+1/nlayers, 1/nlayers)` is brittle at nlayers boundaries (floating-point accumulation can produce `nlayers+1` or `nlayers+2` elements); `np.linspace(0, 1, nlayers+1)` is the idiomatic replacement. The subsequent `sumkh` loop (282-288) is also a straight vectorisable numpy expression.

Not a correctness issue (tests pass), but makes the function harder to read and ~10x slower than it needs to be at large `nlayers`.

### 4. Mixed type signature: `phi_h_fracs: list = None`, then sometimes numpy-array, sometimes list
`layer.py:225, 268, 270, 279, 294`. Inside `lorenz_2_layers`, `phi_h_fracs` flips between Python list (after normalisation at 268-270) and numpy array (line 279 via `np.array`). Downstream `np.array(phi_h_fracs)` at 294 covers both cases but the intent is muddy. Pick one container and stick with it.

### 5. Test coverage gap: no error-case tests in `test_layer.py`
Correctness review notes all five input-validation branches are present (lorenz out of bounds, bad method, bad kh_frac/phih_frac), but no `pytest.raises` tests cover them. Other modules already migrated to `test_errors.py` (per CLAUDE.md tests.md). Layer's 5 validation branches are uncovered.

### 6. Test coverage gap: Langmuir branch of `lorenz_from_flow_fraction` untested
`test_layer.py` tests roundtrips with default EXP method only. The LANG branch of `lorenz_from_flow_fraction` (line 151-155) — the closed-form `B = (y - x) / (x * (1 - y))` — has no test. Add a `test_lorenz_from_flow_fraction_lang` pair.

### 7. Test coverage gap: `B` parameter as pre-computed input
Both `lorenz_2_flow_frac(..., B=...)` and `lorenz_2_layers(..., B=...)` advertise that pre-supplying `B` skips the solver. No test exercises this path — if the short-circuit at line 203 / 255 ever breaks, nothing catches it.

### 8. Test coverage gap: `phi_h_fracs` user-layer path
Three distinct sub-paths in `lorenz_2_layers` — sum>1, sum<1, sum=1 — none tested. Covered in the correctness review manually; should be pinned in test_layer.py.

## Low

### 1. Private helper `bisect_solve` imported from `shared_fns` — bisection limits `(lo, hi)` are per-callsite duplicated
`layer.py:77-80, 158-159`. Both `lorenz2b` and `lorenz_from_flow_fraction` redeclare `hi = 709 (EXP)` / `hi = 25000 (LANG)` / `lo = 0.000001`. Once constants are named (high item #2), this duplication collapses automatically.

### 2. `rtol = 0.0000001` hard-coded twice
`layer.py:98, 176`. Convergence tolerance repeated. Fine to leave, but if ever retuned, two places to change. Extract `_BISECT_RTOL = 1e-7`.

### 3. Inner closures `LorenzErr` / `BErr` take `args` as a positional tuple
`layer.py:83-96, 162-174`. Pattern is `def LorenzErr(args, B): lorenz, method = args`. Works, but reads oddly — `bisect_solve` could take the callable as `partial(LorenzErr, lorenz=lorenz, method=method)` and accept `(B,)` directly. Cosmetic; not worth a refactor unless `shared_fns.bisect_solve` is touched for another reason.

### 4. `__init__.py` is empty (0 lines)
Most modules in this codebase use `from .<module> import *`. `layer/__init__.py` having 0 bytes works (the package is importable) but prevents `from pyrestoolbox.layer import lorenz2b` without going through `layer.layer`. The top-level `pyrestoolbox/__init__.py` lazy-loader compensates, but this file is a lone outlier in its pattern. Align with the other modules.

### 5. Type hints are partial
- `lorenz_2_layers` return type is `np.ndarray` — fine.
- `lorenz_from_flow_fraction`, `lorenz_2_flow_frac` return `float` — correct.
- But `phi_h_fracs: list = None` has no `Optional[List[float]]` and no element type. Minor.

### 6. Dead `np.random.shuffle` path when `user_layers=True`
`layer.py:295-297`. If `shuffle=True` and user supplied `phi_h_fracs`, shuffle is silently skipped. Documented in the docstring, but the `if shuffle: if not user_layers:` structure makes the dead case a two-line no-op. Either raise a warning ("shuffle ignored because phi_h_fracs was supplied") or flatten.
