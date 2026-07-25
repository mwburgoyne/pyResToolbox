# Tech Debt — plyasunov

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/` — `iapws_if97.py` (152 lines), `plyasunov_model.py` (364 lines), `water_properties.py` (47 lines), `__init__.py` (12 lines).

**No [0] placeholder convention exists** in the coefficient arrays: `P_IN_*` tables are 5x7 with `i = 1..5` handled via a 1-based `i_index` argument passed to `_compute_ai`. A separate `a0` is computed from `B12`. So the expected pad-with-dummy-zero pattern (common in legacy Fortran-style IAPWS implementations) is correctly absent — good.

Module is internal, single-purpose, and the core code is tight. Debt is concentrated in three places: duplicate gas-name aliasing across three dicts, a never-consumed `rho_and_kappa` convenience, and zero dedicated test coverage.

## High

### 1. Gas-name alias duplication across three parallel dicts
`plyasunov_model.py:139, 222-233, 236-247`. The n-butane aliases `NC4H10`, `N-C4H10`, `NC4` are maintained in:
- `B12()` elif chain (line 139) — a single tuple `('NC4H10', 'N-C4H10', 'NC4')`
- `P_IN_MAP` dict (lines 229-231) — three explicit keys pointing at the same list
- `MW_GAS` dict (lines 243-245) — three explicit keys with the same value

Adding or removing an alias requires editing three locations. Extract a single source of truth:
```python
_ALIASES = {'NC4H10': ('N-C4H10', 'NC4')}
def _canonical(gas): ...
```
then have `B12`, `A12_inf`, `gas_mw` all normalise first. Eliminates drift risk.

### 2. `rho_and_kappa` is dead — never imported anywhere
`iapws_if97.py:140-152`. Defined to avoid recomputing `_gamma_derivatives` twice, but no consumer exists. `plyasunov_model.py:336-337` calls `rho_w(T,P)` and `kappa_T(T,P)` separately, each going through `_check_region1` + `_gamma_derivatives`. So the module *pays* the duplication `rho_and_kappa` was designed to avoid. Either wire `_V2_inf_cached` through `rho_and_kappa` (fast path) or delete the function.

Bonus: `rho_and_kappa` does **not** call `_check_region1`, while `rho_if97` and `kappa_T_if97` do. If the function were used, it would silently compute out-of-range values. Either validate or delete.

### 3. No test coverage for the plyasunov package at all
No `test_plyasunov.py` exists. The module is exercised indirectly via `test_brine.py` (through `SoreideWhitson`), but:
- `B12()` for each of the 8 gases — untested per-gas
- `A12_inf()` per-gas — untested
- `V_phi`/`V2_inf` LRU cache behaviour — untested
- `gas_mw()` round-trip and error path — untested
- IAPWS-IF97 `rho_if97`/`kappa_T_if97` numerical values against published references — untested
- Alias normalisation (`NC4H10` vs `NC4`) — untested

Add a `test_plyasunov.py` with known-value spot checks: B12 vs paper tables, V2_inf at a couple of reference T/P states (Plyasunov papers contain example tables), and roundtrip `gas_mw('NC4') == gas_mw('NC4H10')`.

## Medium

### 1. Unused re-exports in `__init__.py`
`plyasunov/__init__.py:12` re-exports `V2_inf`, `B12`, `A12_inf`, `MW_GAS` — none of which are consumed outside the package. Only `V_phi` and `gas_mw` are imported by `brine.py`. Dead public surface.

Also covered as an API/UX finding; from a tech-debt angle the issue is that any refactor of `V2_inf`/`B12`/`A12_inf` signatures must worry about unknown external consumers even though there are none.

### 2. `MW_WATER` is defined in `water_properties.py` and re-imported by `plyasunov_model.py`, but `TC_WATER` is imported-and-re-used with no documentation of the dependency
`plyasunov_model.py:33` imports `rho_w, kappa_T, TC_WATER, MW_WATER` from `water_properties`. `OMEGA = 1e-3 / MW_WATER` at line 38 is the only use of `MW_WATER`; `TC_WATER` is used once at line 263. No comment explains why the water module owns these constants. If a reader touches `water_properties` they have to trace back to `plyasunov_model` to know those constants are consumed. Single-line comment would help.

### 3. `_B12_square_well` uses units inconsistently with rest of file
`plyasunov_model.py:74-89`. `sigma` comes in Angstrom, is converted to cm inline as `sigma_cm = sigma * 1e-8`. Nothing else in the file uses Angstrom; all other constants are SI-ish or cm3/mol. Extract `_ANG_TO_CM = 1e-8` with a comment, or convert the table values `_SW_CH3/_SW_CH2/_SW_H2S` to cm at definition time.

### 4. Magic polynomial/exponent literals in `_B12_CO2`
`plyasunov_model.py:58-71`. The hard-coded exponents `-0.5, 0, -3, -6, -10.5` and coefficient list `b = [15.210, 149.72, ...]` are inlined inside the function body. CLAUDE.md's "Named Constants" rule pushes these to module scope:
```python
# Hellmann (2019), Part II Eq. 8
_B12_CO2_COEFFS = (15.210, 149.72, -534.54, -2234.6, -13017.0, -39482.0)
_B12_CO2_EXPONENTS = (0.0, -0.5, -1.0, -3.0, -6.0, -10.5)
```
Also: `b[0] + b[1]/Tstar**0.5 + ...` is a six-term polynomial — the `zip(exponents, coeffs)` loop is simpler and matches `_B12_polynomial`'s style.

### 5. `_compute_ai` exponent table is inlined every call
`plyasunov_model.py:264-266`:
```python
exponents = [-0.5, -1.0, -2.0, -3.0, -4.0, -5.0]
n_map = {1: 9, 2: 8, 3: 7, 4: 6, 5: 6}
```
Allocated on every call (and `A12_inf` calls `_compute_ai` 5 times per evaluation). Move to module-scope constants. Tiny perf gain, but mostly a style / intent signal.

### 6. `V_phi` vs `V2_inf` — one is advertised as an alias of the other with no functional reason for both to exist
`plyasunov_model.py:342-356`. `V_phi(gas, T, P)` body is `return V2_inf(gas, T, P)`. If retained for naming reasons (Plyasunov paper uses `V_phi`, IAPWS convention uses `V2_inf`), collapse to one `__all__` with a comment explaining the duplication. Currently both names are importable from `__init__.py` → two public API points with identical semantics.

### 7. Coefficient-array docstring mentions "Paper 5" / "Paper 6" / "Paper 7" / "Paper 8" without cross-reference
`plyasunov_model.py:154, 179, 212`. Says "Paper 8 (Part IV, Table 2) - full precision - supersedes Paper 5" but the module docstring doesn't define what "Paper 5/6/7/8" are. A reader has to infer that they track Plyasunov's paper sequence. Either reference the full bibliographic citation once and drop the `Paper N` shorthand, or add a legend at module top.

### 8. `iapws_if97.py` module docstring says "no external dependencies beyond numpy" but the file does not import numpy
`iapws_if97.py:4`. The file uses pure-Python `**` arithmetic — no `numpy` import (`import numpy` is nowhere in the file). Not broken, but the docstring is misleading. Either add `import numpy as np` (if planned) or fix the docstring to "no external dependencies".

## Low

### 1. Type hints absent throughout plyasunov_model
`plyasunov_model.py`. No `-> float` annotations on any of `B12`, `A12_inf`, `V2_inf`, `V_phi`, `gas_mw`, `_compute_ai`, `_B12_polynomial`, etc. Brine code imports `V_phi` and uses return as float — type hints would catch a silent tuple-return regression. Adding `(gas: str, T: float, P: float) -> float` is mechanical.

### 2. `_compute_ai` 1-based `i_index` argument is error-prone
`plyasunov_model.py:254-273`. Callers at line 300 do `_compute_ai(p_in[i], T, i + 1)` where `i` is 0-based. The `+1` offset to match paper indexing is a known footgun — easy to forget on next edit. Either keep 0-based `n_map` + expose it that way, or rename the arg to `paper_index` with a clearer comment.

### 3. `P_IN_MAP` and `MW_GAS` gas-name dict keys include both `NC4H10` and `N-C4H10` but `B12()` elif at line 139 also hard-codes `NC4`
Asymmetry: `B12()` accepts `NC4` via the `in (...)` tuple; `P_IN_MAP` has it as a dict key. Once they're unified per High #1, this goes away.

### 4. `iapws_if97.py:_REGION1_IJN` is a list of 3-tuples iterated linearly every call
`iapws_if97.py:68-94`. `_gamma_derivatives` iterates the 34-entry list on every evaluation. `V_phi` is LRU-cached so this rarely matters, but extracting three `numpy` arrays `(I_arr, J_arr, n_arr)` and vectorising the sum would shorten `_gamma_derivatives` to ~5 lines — and match the style used elsewhere in the codebase (gas BNS EOS arrays). Out of scope for pure tech-debt unless someone touches this file.

### 5. `water_properties.py` functions `rho_w`/`kappa_T` are thin wrappers that just cast inputs to float and call the backend
`water_properties.py:22-47`. The `float(T), float(P)` cast is the only value-add (the wrappers exist to let plyasunov substitute backends). Fine, but the docstring doesn't explain the coercion — if caller passes `np.float64` it is downcast to Python `float`, losing no precision but possibly surprising a numpy-heavy caller. One-line note.

### 6. `MW_WATER = 18.015268` lives in `water_properties.py` but `pyrestoolbox/constants/` also defines a water molecular weight
`plyasunov/water_properties.py:18`. Potential duplication with `pyrestoolbox.constants`. Worth a quick audit: if they agree, import from constants; if they disagree, document why (IAPWS vs IUPAC convention values differ at 6th decimal).

### 7. `_check_region1` redundant with `rho_and_kappa` call sites — see High #2
Already covered above; noted here because it's a quality-of-validation tech-debt angle, not just a dead-code angle.

### 8. Plyasunov Part IV 2021 citation missing from References block
`plyasunov_model.py:19-29`. Referenced inline at line 8, not in the References block. Cosmetic / reproducibility. Also flagged in correctness review.
