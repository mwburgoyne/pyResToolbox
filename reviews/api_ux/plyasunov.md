# API/UX — plyasunov

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/` (`__init__.py`, `plyasunov_model.py`, `water_properties.py`, `iapws_if97.py`).

The module is documented as internal. Verified: **not** listed in `pyrestoolbox/__init__.py` `submodules` (line 44-60), so `pyrestoolbox.plyasunov` is only reachable via explicit import path. The brine module consumes it via the private path `from pyrestoolbox.plyasunov import V_phi, gas_mw` and aliases with a leading underscore (`_plyasunov_V_phi`, `_plyasunov_gas_mw`) — good hygiene, no leak to user surface.

No dedicated RST file exists for plyasunov, and none should — it is internal. No issues there.

## Critical
None.

## High

### 1. Public `__init__` re-exports more symbols than brine actually needs — expands the unintended API surface
`plyasunov/__init__.py:12` exports `V_phi, gas_mw, V2_inf, B12, A12_inf, MW_GAS`. But brine only imports `V_phi` and `gas_mw`. `V2_inf`, `B12`, `A12_inf`, `MW_GAS` are re-exported for no current consumer. Since the module is advertised as internal, this expands what an adventurous user could reach via `from pyrestoolbox.plyasunov import A12_inf`. Trim the re-export to the actual brine dependency, or prefix internals with underscore.

### 2. No `__all__` in either `__init__.py` or `plyasunov_model.py`
`plyasunov_model.py` has no `__all__`. Combined with broad `from .plyasunov_model import ...` in `__init__.py`, there is no gate on what is importable. Private helpers like `_B12_polynomial`, `_B12_CO2`, `_B12_square_well`, `_compute_ai` are already `_`-prefixed (good), but future additions without underscore will silently become "public-ish".

## Medium

### 1. `gas` parameter accepts inconsistent aliases across functions
`plyasunov_model.py:139` `B12` accepts `NC4H10`, `N-C4H10`, `NC4` for n-butane. `P_IN_MAP` (line 222) and `MW_GAS` (line 236) duplicate the alias entries. `V_phi`, `V2_inf`, `A12_inf` route through the same map so they inherit the aliases — good. But the docstring of `V_phi` (line 349) lists only `NC4H10`, not the aliases. If the aliases are supported, document them; if they are not promised to remain, remove them.

### 2. `water_properties.MW_WATER` and `TC_WATER` are imported by `plyasunov_model.py` but the constants section of `plyasunov_model.py:35-38` also defines its own `OMEGA = 1e-3 / MW_WATER`
Not incorrect, but the dependency direction (model imports from water_properties) is easy to miss. A quick note in the module docstring would help.

### 3. `B12`, `A12_inf`, `V2_inf`, `V_phi`, `gas_mw` docstrings use `Parameters:` / `Returns:` colon style instead of numpydoc/Sphinx style used elsewhere in the codebase (`Parameters\n----------`)
Inconsistent with the rest of the package (compare `recommend.py:74-92` or `sensitivity.py:118-137`). Since the module is internal, this has no agent-facing impact, but if any of these ever migrate to the public surface the docstrings need reformatting.

### 4. `V_phi` is documented as "Alias for V2_inf" — unclear why both exist
`plyasunov_model.py:342`. `V_phi` just calls `V2_inf(gas, T, P)`. If the distinction is nominal (symbol convention used in different papers) it is worth stating that explicitly rather than calling one an alias of the other.

### 5. `gas_mw` error message does not list supported gases
`plyasunov_model.py:363`. Raises `ValueError(f"Unknown gas: {gas}")` — whereas `A12_inf` on line 293 raises `ValueError(f"Unknown gas: {gas}. Available: {list(P_IN_MAP.keys())}")`. Inconsistent — `gas_mw` should show the valid set too.

### 6. `B12` error message missing supported-gas list
`plyasunov_model.py:144`. Same issue as above — message is `f"Unknown gas: {gas}"` with no list.

## Low

### 1. Unit string capitalisation is inconsistent in docstrings
"cm3/mol" vs "MPa" vs "kg/m3" — fine — but "Kelvin" spelled out once in `water_properties.py:11` and abbreviated "K" everywhere else. Purely cosmetic.

### 2. `B12` returns `cm3/mol` but the convention in brine module is SI — callers handle conversion externally
Not a problem per se, but worth a one-line note that `V_phi` returns cm3/mol (Plyasunov convention) since brine code has to divide/multiply for kg/m3 density corrections.

### 3. `rho_w` and `kappa_T` (in `water_properties.py`) wrap `rho_if97` and `kappa_T_if97` by name only
These are trivial passthroughs that cast to `float(T), float(P)`. Fine; the wrappers exist to give the public surface neutral naming (not tied to IF97). Document that in the file docstring — currently only the backend is mentioned.

### 4. Module docstring claims "No external dependencies" for water_properties
`water_properties.py:8` says "no external dependencies" — but `iapws_if97.py` uses `numpy`. The intent was probably "no third-party IF97 library like `iapws`" — reword to avoid confusion.

### 5. No bibliographic entry for Plyasunov Part IV
Noted also in correctness review (item 4). Part IV 2021 (H2/N2/CH4) is referenced inline (`plyasunov_model.py:7-8`) but missing from the References block (lines 19-29). Since this module is internal, impact is purely for developer agents debugging the implementation.
