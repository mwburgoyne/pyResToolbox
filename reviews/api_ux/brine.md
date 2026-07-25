# API/UX — brine

Scope: `pyrestoolbox/brine/` (brine.py, _lib_vle_engine.py, _lib_salting_library.py) + externally-facing surface of `pyrestoolbox/plyasunov/`. Correctness items out of scope (see `reviews/correctness/brine.md`). Prior 2026-04-12 review Tier 4 items (parameter alias standardization, VLE damping, V_phi caching) already done — this review verifies consistency only.

---

## Critical

- **`SoreideWhitson` hides `framework` and `salinity_method` from end users** — `brine.py:1738-1741`. The underlying VLE engine has two dispatch parameters that materially change numerical results: `framework ∈ {'proposed', 'sw_original', 'dropin'}` and `salinity_method ∈ {'gamma_phi', 'explicit', 'embedded'}`. The top-level `SoreideWhitson.__init__` hardcodes `framework='proposed'`, `salinity_method='gamma_phi'` with no way to override. The correctness review notes `framework='proposed' + salinity_method='embedded'` silently falls back to `gamma_phi`; here the issue is different — users cannot even select a framework. The class docstring (brine.py:1406-1477) and `brine.rst:310-489` make no mention of these knobs. A user wanting to reproduce published S&W 1992 results, or to compare frameworks, must bypass the public API entirely. Either (a) expose both as constructor kwargs with defaults, or (b) explicitly document "proposed framework only" as a design decision.

---

## High

- **`brine_props` return shape is non-self-describing and asymmetric** — `brine.py:384`. Returns `(bw, lden, visw, [cwu, cws], rsw)` — a 5-tuple where index [3] is a 2-element list. No named-tuple, dataclass, or dict wrapper. Users must know positional order. `CO2_Brine_Mixture` and `SoreideWhitson` return class instances with attributes — the three models' output shapes are inconsistent. Consider a `BrineProps` NamedTuple or dataclass, or at minimum document the return-shape rationale (legacy API compat).

- **`ppm` vs `wt` alias auto-conversion can silently override user intent** — `brine.py:136-137, 482-483, 1486-1487`. `brine_props`, `CO2_Brine_Mixture`, and `SoreideWhitson` all do `if wt is None: wt = ppm / 10000 if ppm is not None else 0` (or the inverse). But if a user passes both `wt=5` and `ppm=50000`, the explicit kwarg wins silently — no check for inconsistency, no warning. Worse: passing both `wt=5` and `ppm=10000` (which disagree) will pick whichever comes first in the resolution chain with no feedback. Recommend: raise `TypeError` when both are specified.

- **Error message for invalid gas fraction sum is misleading** — `brine.py:1681-1683`. When `y_hc < 0` (non-HC fractions sum > 1), the error reads `f"Non-HC gas fractions sum to {1 - y_hc:.4f}, exceeding 1.0"`. But `1 - y_hc` is not the sum — it's `1 - (1 - sum)` = sum + (negative_term). Actually when `y_hc = 1 - y_CO2 - y_H2S - y_N2 - y_H2` is negative, `sum = 1 - y_hc` is > 1, so the arithmetic is correct but the construction is confusing. Fine on inspection but the value reported doesn't match intuition — use `y_CO2 + y_H2S + y_N2 + y_H2` directly. Also, constructor at line 1498-1500 has its own separate check raising a different message for the same condition — pick one.

- **`SoreideWhitson.y` attribute is misleading** — `brine.py:1745`, `brine.rst:394-395`. `.y` is documented as "Gas phase compositions (dry basis, normalized)", but it is set to `dict(self.gas_comp)` — the input-side normalized gas composition, NOT the computed vapor-phase composition from the flash (that's `.y_H2O` + `(1 - y_H2O) * self.gas_comp` implicitly). Compare to `CO2_Brine_Mixture.y = np.array([yCO2, yH2O])` which actually is the vapor composition including water. Users will be confused — `.y` does not include the water fraction. Either rename (`.y_dry`, `.gas_composition_dry`) or return the full wet-gas composition dict.

- **`.Rs` output type differs across models** — `CO2_Brine_Mixture.Rs` is a float, `SoreideWhitson.Rs` is a dict (`brine.rst:186-188` vs `:415-417`). This is a real semantic difference (multi-gas needs per-gas breakdown), but it means any user code written for `CO2_Brine_Mixture` cannot be trivially ported to `SoreideWhitson`. Consider making `CO2_Brine_Mixture.Rs` a dict too (`{'CO2': ...}`) for API uniformity, keeping float as `Rs_total`.

---

## Medium

- **Default `metric` flag differs across the three models** — `brine_props` defaults `metric=False` (oilfield), but `CO2_Brine_Mixture` and `SoreideWhitson` default `metric=True`. This is called out in `brine.rst:17-21` but it's still a landmine: a user switching from one model to another will get silent unit shifts. CLAUDE.md (Tier 1) says "All APIs default to oilfield units." — two of three brine classes violate this. Either align all three to `metric=False`, or accept the inconsistency and surface it more prominently (maybe require explicit `metric=` to avoid ambiguity for the CO2/SW classes).

- **Component naming convention is inconsistent** — the VLE engine uses VLE-specific keys (`'CH4'`, `'C2H6'`, `'nC4H10'`, `'H2S'` — all uppercase with chemical-formula syntax), but `SoreideWhitson` constructor uses argument names (`y_CO2`, `y_H2S`, `y_N2`, `y_H2`) — capital letters in kwargs. `.Rs` / `.x` / `.gas_comp` dicts then key back to VLE-style (`'CH4'`, `'nC4H10'`). No alternate spellings supported: `'co2'`, `'CarbonDioxide'`, `'ch4'`, `'methane'` all fail. For AI agents this isn't discoverable without reading source. Document valid component keys in the class docstring; consider adding case-insensitive lookup or at least an error message that enumerates valid keys when an unknown gas appears (currently `brine.py:1681` only handles the sum condition, not unknown gases — `KIJ_AQ_PROPOSED` at `_lib_vle_engine.py:602` does list valid keys but that path isn't hit from `SoreideWhitson.__init__`).

- **Many internal VLE engine classes and functions are accessible but undocumented** — `SWBinaryVLE`, `SWMultiComponentFlash`, `H2WaterVLE`, `calc_gas_brine_equilibrium`, `kij_aq_*` correlation fns are imported via `from pyrestoolbox.brine._lib_vle_engine import ...` (underscore prefix marks the file as private, good) but there's no `__all__` restriction in `_lib_vle_engine.py`. Users who want direct VLE access have to import from a leading-underscore module — which Python convention says is private. Either expose a curated public VLE API (e.g. `brine.vle_engine.SWBinaryVLE`) or explicitly document that `_lib_vle_engine` is internal and offer a narrow facade.

- **`SoreideWhitson` docstring example is inconsistent with code** — `brine.py:1456`, docstring says `mix.Rs  # Returns per-gas Rs dict, e.g. {'CO2': 15.2}` without specifying units; `brine.rst:444-445` shows `{'CO2': 140.90858294709142}` for a pure-CO2 5000-psia case. 140 scf/stb is reasonable. The discrepancy between docstring example (15.2) and RST example (140.9) for comparable conditions will confuse users — align them.

- **`make_pvtw_table` is a deprecation wrapper but the docstring doesn't say "deprecated" clearly** — `brine.py:1349-1352`. Docstring says "Deprecated: Use simtools.make_pvtw_table() instead." but doesn't emit `DeprecationWarning`. RST at :221-224 mentions it's a wrapper without marking `.. deprecated::`. If this is truly deprecated, add the warning or drop "Deprecated" word from the message.

- **`plyasunov` submodule leaks internals** — `pyrestoolbox/plyasunov/__init__.py` exports `V_phi`, `gas_mw`, `V2_inf`, `B12`, `A12_inf`, `MW_GAS`. Tier 2 rule says this submodule is "internal — not exposed via `__init__`" of the main package — correct, it's not in `pyrestoolbox.__init__`. But the submodule is still importable as `pyrestoolbox.plyasunov`, and its `__init__.py` exposes six names without a `__all__`. If it is internal, consider prefixing exports with underscore or documenting the intended surface. `V2_inf`, `B12`, `A12_inf` look like low-level correlation internals that external users shouldn't touch.

- **`SoreideWhitson` constructor accepts per-gas mole fractions but not a composition dict** — API is `y_CO2=0, y_H2S=0, y_N2=0, y_H2=0`. For an AI agent with a composition dict `{'CO2': 0.1, 'N2': 0.03, 'CH4': 0.87}`, this requires manual unpacking and silently drops HC components (C2+) because they're not exposable via kwargs (user relies on SG-based HC splitter). Consider accepting an optional `gas_composition: dict` kwarg that overrides both the y_* kwargs and the SG-derived HC split.

---

## Low

- **`CO2_Brine_Mixture` docstring has a stale `xSalt` note** — `brine.py:450`: "Sum of Mole fraction of Na + Cl species in brine (Double the NaCl species mole fraction)". This is correct but the parenthetical can confuse — users may expect `.xSalt` to be xNaCl (half of what it actually is). `brine.rst:169-170` says "Mole fraction of NaCl in brine" which contradicts the code comment. Align the two.

- **`CO2_Brine_Mixture` has 30+ attributes set in `__init__` as None placeholders** — `brine.py:491-529`. Users doing `dir(mix)` get a lot of noise (`.fugPi`, `.K`, `.gamma_prime`, `.pRT`, `.A`, `.Bprime`, `.aMix`, `.aij`, `.kij`, `.bMix`, `.b`, `.vBar`, `.CO2_sat`, `.repeat`, `.scaled`, `.low_temp`, `.P0`, `.pRT0`, `.Rs_STD`, `.ppm_sat`, `.tKel`, `.degC`, `.pBar`, `.molaL`, `.MolarVol`, `.GASZ`, `.MwBrine`, `.MwGas`). These are all iteration state, not user-facing output. Prefix with `_` or move into a `_state` dict.

- **`CO2_Brine_Mixture` docstring inputs section uses `pres`/`temp`/`ppm` only** — `brine.py:441-445` and RST `:138-152` don't mention the `p`/`degf`/`wt` aliases, even though `__init__` accepts them. Same for `SoreideWhitson` (`brine.rst:350-358`). Users reading docs won't know aliases exist. Add a note in each inputs table.

- **`brine_props` docstring documents the tuple format inline** — `brine.py:123-124`: "Returns Tuple of (Bw (rb/stb | rm3/sm3), Density (sg), viscosity (cP), Compressibility (1/psi | 1/bar), Rw GOR (scf/stb | sm3/sm3))". Item [3] is a list, not a scalar. The docstring collapses it to "Compressibility (1/psi | 1/bar)" without mentioning it's a 2-element list. RST at :93-95 is correct. Sync the docstring.

- **`SoreideWhitson.Cf_sat = None` when `cw_sat=False`** — `brine.py:1561`. Returning `None` silently for an unset calculated value is a common pattern but annoying — callers that unconditionally do arithmetic on `Cf_sat` will get TypeError rather than a helpful "call with cw_sat=True" message. Consider a lazy property or a sentinel object with a descriptive repr.

- **`CO2_Brine_Mixture.bVisblty` vs `SoreideWhitson.bVisblty`** — both exist, same meaning, typo ("viscosibility" is fine as a coined term but the truncation "bVisblty" is non-obvious). Probably too late to rename, but note it. Neither is exported in `__all__` as a canonical accessor.

- **`brine/__init__.py` is just `from .brine import *`** — relies on `__all__` in `brine.py:36-38` which lists `['brine_props', 'CO2_Brine_Mixture', 'SoreideWhitson', 'make_pvtw_table']`. Fine, but the `__init__.py` has no module-level docstring. Helps agentic discoverability if `brine/__init__.py` had a one-line summary and explicit re-exports (`from .brine import brine_props, CO2_Brine_Mixture, ...`) rather than wildcard import.

- **The docstring in `SWMultiComponentFlash.calc_equilibrium`** (`_lib_vle_engine.py:2249-2251`) documents `salinity_for_sechenov` as "(Deprecated)" but it's still accepted. Since this is an internal API, either remove the deprecated param or emit `DeprecationWarning`.

- **`calc_gas_brine_equilibrium` default `method='henry'`** (`_lib_vle_engine.py:2369`) while `SoreideWhitson` hardcodes `method='flash'` (`brine.py:1738`). Unsuspecting users calling `calc_gas_brine_equilibrium` directly get Henry's law binary decomposition (valid only for dilute) instead of the more accurate flash. Recommend: default `method='flash'` at the engine level; the existing flash-with-adaptive-damping is robust enough.

---

## Summary

Critical: 1 — High: 5 — Medium: 7 — Low: 9

Top 3:
1. **`framework` and `salinity_method` hidden from `SoreideWhitson` users** — the most impactful VLE knobs are literally not exposable without going into `_lib_vle_engine`.
2. **`.y` attribute semantics wrong** — `SoreideWhitson.y` returns the input dry-gas composition, not the flash vapor composition, despite being documented as the latter and despite `CO2_Brine_Mixture.y` meaning the vapor phase.
3. **`.Rs` type asymmetry (float vs dict) across the three models** — makes model-swapping code non-portable; worth a uniform dict-of-components interface.
