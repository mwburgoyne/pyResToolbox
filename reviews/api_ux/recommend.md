# API/UX — recommend

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/recommend/recommend.py`, `/home/mark/projects/pyResToolbox/pyrestoolbox/recommend/__init__.py`, `/home/mark/projects/pyResToolbox/pyrestoolbox/docs/recommend.rst`.

Module is tidy: 4 functions + 1 dataclass, clean `__all__`, numpydoc docstrings, and a well-aligned RST. The main issues are two dead parameters (`sg`, `well_type`) and soft-hardcoded alternatives lists that duplicate information available from the Enum classes.

## Critical
None.

## High

### 1. `sg` parameter in `recommend_gas_methods` / `recommend_methods` is never used
`recommend.py:71` and `recommend.py:237`. The function accepts `sg: float = 0.65` and threads it through `recommend_methods` (line 268), but no branch in the decision logic references it. Users will reasonably expect that specifying a heavy gas (`sg=0.9`) changes the recommendation (e.g., DAK vs BNS — DAK has looser validity at high sg), but it does nothing. Either wire it into the decision tree or remove the parameter. Documented in both the function docstring and the RST (`recommend.rst:32`) with no caveat — this is user-facing misinformation.

### 2. `well_type` parameter in `recommend_vlp_method` / `recommend_methods` is never used
`recommend.py:199` and `recommend.py:241`. Same pattern as `sg` — declared, documented, but no branch consumes it. Flagged in correctness review too. This one is arguably more misleading because a user would reasonably expect gas vs oil to change VLP recommendations (HB was specifically developed for oil/gas mixtures, Gray for gas condensate — the recommendation does change in practice, just not in this code).

## Medium

### 1. `alternatives` lists are hand-curated — drift risk against `classes.py` Enums
`recommend.py:99, 116, 121, 129, 134, 142, 147, 181, 186, 192, 223, 231`. All hard-coded strings (`'DAK'`, `'PMC'`, `'WG'`, etc.). If a new method lands in `classes.z_method` or `vlp_method`, the recommend engine will not surface it without a matching manual update here. Consider pulling candidate lists from the Enum classes or cross-validating in a test.

### 2. No input validation — negative fractions or inert sum > 1 silently accepted
`recommend.py:71-93`. Other gas-module entry points call `validate_pe_inputs`. Here, `recommend_gas_methods(co2=0.8, h2s=0.8)` returns a "High inert content (160%)" recommendation with a sensible-looking rationale string. Covered in correctness review, but worth re-stating as a UX issue: the user is not told that their composition is nonsense.

### 3. `recommend_oil_methods` gives the same top recommendation for every API
`recommend.py:154-195`. Across `api < 10`, `10 <= api <= 50`, `api > 50`, all branches return `VELAR` / `VELAR` / `MCAIN`. Only the `rationale` string changes. If this is truly the design (VELAR is universally preferred), the function does not need the if/else scaffolding and the presence of that scaffolding implies tiered recommendations that do not exist. Simplify or actually differentiate.

### 4. Thresholds (55% inerts, 10% CO2/H2S, 30° deviation) are magic numbers
`recommend.py:111, 124, 216`. Not wrong, but the CLAUDE.md coding guide ("Named Constants") pushes for `_INERT_THRESHOLD_BNS = 0.55` etc. with paper citations. Low-friction fix.

### 5. `MethodRecommendation.alternatives` is an ordered list with no documented meaning
Is the first element the "next-best choice"? Alphabetical? Priority-ranked? The docstring says "Other viable methods" but downstream users will be unsure whether to pick `alternatives[0]` or any. Document the ordering convention (the code appears to use priority order, but that is not stated).

### 6. Master `recommend_methods` silently drops oil recommendations when `api=None` but returns a dict with different shapes
`recommend.py:270-273`. The return dict has 3 keys when `api=None` and 6 keys when `api` is set. That is unusual for an API — a caller writing `recs['pbmethod']` gets `KeyError` when they forgot to pass `api`. A safer pattern would be always-6 keys, with missing values as `None`. Documented in RST (`recommend.rst:247-254`) as "only when api is provided", so it is not a silent surprise — but dict shape shifts are a UX footgun.

## Low

### 1. `rationale` strings format percentages with `{h2:.1%}` even when `h2 = 0`
Does not currently fire because the H2 branch gates on `h2 > 0`, but if the inert branch text ever includes an individual percentage for a zero-value component, it would say "0.0%". Minor.

### 2. `MethodRecommendation` dataclass has no `__repr__` customization
Default dataclass repr will print the full rationale string on one line. Not awful, but for a print-heavy use case (`print(recommend_methods(...))`) a multi-line formatter would help. Out of scope if simplicity is preferred.

### 3. Docstring in `MethodRecommendation` lists `mandatory` but never explains when it is `True`
`recommend.py:48-68`. Would be clearer: "`True` when no alternative methods exist (currently only the H2-present case)."

### 4. `__all__` and `__init__.py` are clean
`recommend/__init__.py` is a single-line `from .recommend import *`. `recommend.py:37-41` declares `__all__` explicitly, re-exporting `MethodRecommendation` as a public class. Good.

### 5. RST uses `'zmethod'` / `'cmethod'` as dict keys — matches the Enum class naming convention
No inconsistency between the recommend RST and the gas RST. Good.

### 6. Examples in RST are concise and doctest-friendly
`recommend.rst:63-80, 120-132, 171-185, 256-268`. All examples are executable and value-stable. Good.
