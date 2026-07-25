# Tech Debt — brine + plyasunov

Scope: `pyrestoolbox/brine/` and `pyrestoolbox/plyasunov/`. Tier-4 items already shipped in v3.1.4 (adaptive VLE damping, V_phi/kij caching, parameter aliases) are not re-flagged. Correctness items from `reviews/correctness/brine.md` and API/UX items from `reviews/api_ux/brine.md` are not duplicated.

## High

- **File size / monolith: `_lib_vle_engine.py` is 2726 lines.** Single file holds ComponentProperties database, 3 kij dispatch tables (`KIJ_AQ_PROPOSED`, `KIJ_AQ_SW_ORIGINAL`, `KIJ_AQ_DROPIN`), ~40 kij functions, alpha functions, Sechenov routing, PR-EOS flash scaffold, Rachford-Rice solver, `SWMultiComponentFlash` class, `flash_tp`, `calc_equilibrium`, backward-compat aliases, and a 120-line `__main__` demo. No sub-module boundaries. Navigation, diff review, and targeted test loading all suffer. Split into: `_vle_components.py` (ComponentProperties + constants), `_vle_bip.py` (kij dispatch + all per-gas functions), `_vle_alpha.py` (alpha functions), `_vle_flash.py` (engine + Rachford-Rice + flash_tp), `_vle_api.py` (`calc_equilibrium`).

- **File size: `brine.py` is 1913 lines.** Mixes three independent model classes (`brine_props` function, `CO2_Brine_Mixture`, `SoreideWhitson`) plus `make_pvtw_table` wrapper. `SoreideWhitson` alone is ~500 lines. Splitting into `_brine_ch4.py`, `_brine_co2.py`, `_brine_soreide_whitson.py` would mirror the oil-module refactor already done (Tier 3 roadmap) and reduce import blast radius.

- **Python/Rust parallel flash paths with no parity tests.** `flash_tp` (2106-2116) and `calc_equilibrium` (2270-2305) both have try-Rust / fall-back-to-Python branches; `CO2_Brine_Mixture.co2BrineSolubility` similarly dispatches `_native.co2_brine_solubility_rust`. The Python-fallback path is silently taken whenever the Rust module is unavailable, but `tests/test_brine.py` exercises only the default (Rust-preferred) path. There is no test that forces Python and Rust separately and compares. Correctness review already flagged the `proposed`-framework ks divergence; the deeper debt is that no CI gate catches *any* future Python↔Rust drift. Add a `test_brine_rust_parity.py` that sweeps a grid and asserts `|rust - python| < tol` per output field.

- **Spivey salt-correction coefficient lookup duplicated between `brine_props` and `brine_props_co2`.** `brine.py` defines Spivey coefficient arrays (`_RHOW_T70_ARR`, `_EWT_ARR`, `_FWT_ARR`, `_DM2T_ARR`..`_FM12T_ARR`) once, but the `_Eq41` evaluation logic (lines 147-166 in `brine_props`, duplicated again in `brine_props_co2` around lines 1190-1201) plus the density chain (rhowtp · salt_ratio · Garcia) is re-implemented. Extract to a private helper `_spivey_salt_density(T_K, P_MPa, molality) -> (rho_spivey_tp, rhow_t70, rhow_tp, salt_ratio)` and call it from both entry points. Same applies to the `Ebtm/Fbtm` pressure correction block.

- **`_calc_properties` calls `brine_props` four times for the same (p, T).** In `SoreideWhitson._calc_properties` (lines 1756-1780) `brine_props(p, degf, wt, ...)` is invoked for base density, then again for standard-condition density, then twice more inside the viscosity/compressibility path. Each call re-runs IAPWS-IF97 + Spivey + CH4 solubility — expensive and wasted. Cache the base-brine result once; reuse for all downstream corrections.

## Medium

- **15-gas `kij_aq_*` wrapper bloat.** Lines 240-297 of `_lib_vle_engine.py` hand-write `kij_aq_c2h6`, `kij_aq_c3h8`, `kij_aq_nc4h10`, ..., `kij_aq_nc10h22`, each a three-line trampoline to `kij_aq_hydrocarbon(T, m, COMP[...])`. The `_kij_aq_hc_proposed(gas)` closure factory directly above already shows this can be a loop. Replace the 15 explicit wrappers with `{gas: partial(kij_aq_hydrocarbon, comp=COMP[gas]) for gas in HC_GASES}` and populate the dispatch tables from the dict. Cuts ~60 lines and removes a rename-drift hazard.

- **Duplicated physical constants across three files with different precisions.**
  - `brine.py`: `MWSAL = 58.4428`, `MWWAT = 18.01528`
  - `_lib_vle_engine.py`: `MW_NACL = 58.44`, `MW_H2O = 18.015`
  - `plyasunov/water_properties.py`: `MW_WATER = 18.015268`, plus unit `MW_GAS` dict in brine.py

  Three different MW_water values (18.015 / 18.01528 / 18.015268) and two MW_NaCl values in the same package. Pick one source (`pyrestoolbox.constants`), re-export, and delete the duplicates. Any numerical sensitivity difference should be caught by existing frozen baselines.

- **Unit-conversion helpers redefined in engine.** `_lib_vle_engine.py` defines `fahrenheit_to_kelvin`, `celsius_to_kelvin`, `kelvin_to_celsius`, `kelvin_to_fahrenheit`, `psia_to_pascal`, `bar_to_pascal`, `pascal_to_bar`, `pascal_to_psia` (lines 1147-1306). `pyrestoolbox.constants` / shared utilities already cover these. Dead duplication — delete and import.

- **Three kij dispatch tables with overlapping entries.** `KIJ_AQ_PROPOSED`, `KIJ_AQ_SW_ORIGINAL`, `KIJ_AQ_DROPIN` (lines 455-509) each map all 15 gases but many entries are identical (e.g. CH4 MC-3 form appears verbatim in proposed and dropin). A single registry of `(framework, gas) -> fn` with shared defaults plus framework-specific overrides would make "what's different about framework X?" answerable by inspection instead of diffing three 20-row dicts.

- **Commented-out Ezrokhi density block (~50 lines) in `brine.py` lines 594-642.** Dead code with no TODO or dated comment. Either resurrect behind a feature flag or delete.

- **`_lib_vle_engine.py` backward-compat aliases (lines 2508-2602).** `BinaryVLE`, `H2WaterVLE`, `kij_aq_rational`, `kij_aq_linear`, `kij_aq_chabab_2023`, `kij_aq_lopez_lazaro_2019`, `sechenov_TO`. None appear in the public `brine.__init__.py` surface. If they are purely research aliases from the port, move them to a `_legacy.py` sibling or delete with a deprecation comment in the changelist.

- **`__main__` demo blocks ship inside the wheel.** `_lib_vle_engine.py` lines 2607-2726 (120 lines of example script) and `_lib_salting_library.py` lines 544-693 (150 lines). Both are `if __name__ == "__main__":` guarded so they do not execute on import, but they bloat the installed package (~270 lines of dead bytecode compiled into `__pycache__`). Move to `pyrestoolbox/docs/notebooks/` or an `examples/` subdir outside the package.

- **`CO2_Brine_Mixture.brine_props_co2` is a method that never reads `self`.** Should be a module-level private function `_brine_props_co2(...)`, or `@staticmethod` at minimum. Same smell invites accidental hidden state coupling.

- **30+ `self.x = None` placeholder attributes in `CO2_Brine_Mixture.__init__` (lines 491-529).** Dataclass-style init that then gets populated later by `_calc_properties`. Makes the public attribute surface unclear (user reading `__init__` sees `None`s, not the actual schema). Either use `@dataclass` with explicit types or document which attributes are guaranteed populated post-construct.

- **Plyasunov `P_IN_MAP` / `MW_GAS` duplicate-key aliases.** `NC4H10`, `N-C4H10`, `NC4` all map to the same array in `plyasunov_model.py`. Harmless but dead-weight — a single name-normalizer (strip hyphens, uppercase) called at lookup entry would let the maps hold one canonical key per gas.

- **`_compute_ai` hard-coded `n_map = {1: 9, 2: 8, 3: 7, 4: 6}` magic numbers** (`plyasunov_model.py`). These are the summation upper limits per row of the P_IN matrix. Comment references Plyasunov Part IV Eq 32/34; extract as `_PIN_ROW_LIMITS = (9, 8, 7, 6)` with citation or derive from matrix shape.

## Low

- **1-based indexing via `[0]`-padded arrays undocumented.** `_RHOW_T70_ARR`, all `_EWT/FWT` arrays, `_DM*_ARR`, `_FM*_ARR`, `_DUAN_A/B/C/U/LAMBDA/ETA`, `_IAPWS_VAP_A`, `_MAODUAN_D` all begin with `[0, ...]` so that McCain/Spivey Eq-4.1 index `a[1]..a[n]` works. No module-level comment explains this. One line: `# Leading 0 preserves 1-based indexing from McCain 2009 Table 4-2.` at the first occurrence.

- **`_import_dependencies()` no-op stub.** Empty placeholder from an earlier refactor pass — delete.

- **`_last_vlle_warning` flag is set but never read/reset.** Walks through the module unused — either wire into `calc_equilibrium` diagnostics or remove.

- **`calc_equilibrium` still accepts deprecated `salinity_for_sechenov` kwarg** (lines 2249-2266) with a fallthrough path. Schedule removal after one release; add `DeprecationWarning` now.

- **`ks_akinfiev_h2s_from_tables` has NaN at T=298.15 K** in the interpolation table. Documented in correctness review, noted here as the table itself lives in `_lib_salting_library.py` — a clean fix is a one-line table correction.

- **Gas-name normalization no-op: `.upper().replace('N','N')`** at `_lib_salting_library.py` line 350. Copy-paste of an earlier `replace('-','')`. Cosmetic but broken-looking.

- **`water_properties.py` is a 47-line trivial wrapper** around `iapws_if97.py`, adding only `MW_WATER`/`TC_WATER` constants. Merge into `iapws_if97.py` or move the two constants into `plyasunov/__init__.py`.

- **`plyasunov` submodule has no public `__init__.py` re-exports.** Consumers inside `brine/` reach in with `from pyrestoolbox.plyasunov.plyasunov_model import V_phi`. Minor but a thin `__init__` would formalize the internal API.

- **Spycher-Pruess 99–109 C blending weight is an inline magic expression** in `CO2_Brine_Mixture.co2BrineSolubility` rather than a named helper `_spycher_blend_weight(T_C)`. Cited as Spycher-Pruess 2010 Table 1 — lift to named constant block with citation.

- **`brine_props` `_Eq41` helper** is defined at module scope but only ever called from `brine.py`; fine, but has no docstring referencing Spivey Eq 4.1. One-line cite.

- **Tests cover scalar happy paths only.** `tests/test_brine.py` has no parametrized sweeps across framework × salinity_method (even for valid combos), no Rust-off path, no end-to-end "xCO2 high + high-T" round trip. Correctness review lists the gaps — repeated here as a tech-debt item because the current coverage under-protects future refactors of the items above.
