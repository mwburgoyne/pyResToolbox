# API/UX — simtools

## Critical

- **`make_vfpprod` / `make_vfpinj` return BOTH a dict AND may write a file via `export=True`.** Docs say the function "returns a dict" but side-effecting file writes to CWD are buried under `export=False`. No `filename=` path directory control — the file is written relative to CWD. For agent use this is a footgun: agents can't choose an absolute output path without passing a full path in `filename`, and the docstring does not state that `filename` accepts a full path. Recommend documenting that `filename` may be absolute, and recommend returning the written path in the dict (currently `"eclipse_string"` is returned but no `"written_path"` key).

- **VFP `eclipse_string` vs file-write inconsistency.** The `make_vfpprod` / `make_vfpinj` dict key is `"eclipse_string"` (plain str — list-of-lines style is NOT exposed), whereas `make_bot_og` / `make_pvtw_table` / `rel_perm_table` / `influence_tables` silently write `.INC` files when `export=True` but do **not** return the keyword string at all. Inconsistent: one family returns the string and optionally writes, the other only writes. Agents reading an in-memory keyword need to re-load it from disk for BOT/PVTW/SWOF/AQUTAB.

## High

- **`zip_check_sim_deck` uses `input()` prompts.** It prompts interactively if `files2scrape=[]` (the default) and also if files are missing (`"Continue to zip even with missing files?"`). Calling this in an agent context without TTY deadlocks. Needs a `non_interactive=True` flag, or raise cleanly when stdin isn't a tty. Docstring does not warn. Also has `print(...)` calls gated on `console_summary` but prompts are unconditional.

- **`ix_extract_problem_cells` also uses `input()`** when >1 PRT file exists (line 127). Same agent-hostile pattern. Add explicit filename requirement or `select_first=True` default.

- **`make_vfpprod` / `make_vfpinj` swallow exceptions silently** (line 1377, 1641) replacing failed BHP with `1e-6` and `print(...)` warning. A Mscf/d table silently padded with 1e-6 BHP produces a deck that simulates but gives garbage. `n_failed` is in the return dict but no `warnings.warn`. Should emit a real warning (not `print`) and the dict should carry failure indices so agents can filter.

- **`make_bot_og` delegates to `oil.make_bot_og()` (and `make_pvtw_table` to `brine.make_pvtw_table()`), but these are re-exports, not wrappers.** The signatures must be kept in lock-step manually. Any drift (e.g. the `vis_frac` kwarg added to simtools but missing elsewhere) silently breaks users following the old call site. Worth a single source of truth via `__all__` re-export at the API boundary, not two independent function defs.

- **Parameter inconsistency: `krtable` vs `flo_type` casing.** `rel_perm_table` accepts `krtable='SWOF'` (upper). `make_vfpprod` accepts `well_type='gas'` (lower) but `flo_type='WAT'` (upper). `make_vfpinj`'s `flo_type='WAT'`. Mixing cases within the same module is inconsistent and the code forces `.upper()` / `.lower()` asymmetrically — documented casing in the docstring is the canonical one users will copy. Pick one.

- **`make_vfpprod` SGWFN table name typo in Enum vs docstring.** Docstring line 414: `krtable='SGFN'` but code only accepts `kr_table.SGWFN` (`'SGWFN'`). Docstring line 448 in RST has the same `SGFN` typo. Users copying from docs will hit `ValueError`.

## Medium

- **`rel_perm_table` has 28 parameters, flat.** Corey needs `no, nw, ng`; LET needs 9; Jerauld needs 6. Users must read the body to know which subset applies to their family. The `fname == "COR" / "LET" / "JER"` validation catches this at runtime, but there is no grouped dict input (e.g. `corey_params={'no':2.5,'nw':1.5}`). Not a blocker but noisy for LLM callers. Default `=1` for every parameter silently "works" on wrong families.

- **`fit_rel_perm_best` `all_results` key is a dict of individual fit results**, but this duplicates ~3x the per-family data (residuals arrays, kr_pred arrays). Undocumented memory cost; consumers who only want `best['params']` still ship the other two. Fine as-is, but needs a docstring note.

- **`make_vfpprod` / `make_vfpinj` `vlpmethod='WG'` string default** is inconsistent with the `classes.vlp_method` Enum that everything else in the module uses for method dispatch. All other rel-perm methods document "str or Enum", VFP only documents "str".

- **`make_bot_og` accepts 7 separate method Enums** (`comethod, zmethod, rsmethod, cmethod, denomethod, bomethod, pbmethod`). Each with cryptic string codes (`EXPLT`, `PMC`, `VELAR`, `SWMH`, `MCAIN`). The docstring does **not** list which options are valid for each method. RST also omits this. An agent cannot discover valid strings without reading `classes.py`. `make_bot_og`'s docstring needs the option list inline, or at minimum `See classes.py for options`.

- **`make_vfpprod` oil well default `api=45`** (condensate) contradicts expectation for oil wells where 30–35 is more typical. Single default cannot suit both well types — likely needs per-well-type defaults, or split into `cond_api` / `oil_api`.

- **`rr_solver` return is a 5-tuple (`N_it, yi, xi, V, L`)** — no named keys. Other flash solvers in brine `_lib_vle_engine` return structured dicts. Inconsistent across the package.

- **`rr_solver` docstring doesn't say it raises ValueError** on empty arrays / mismatched lengths / `zi_sum <= 0`. RST example doesn't document any error semantics either.

- **`influence_tables` hits `print("Calculating ReD = ...")` on Python fallback** (line 768). No silent option. Rust path is silent. User sees stdout spam only when Rust isn't available.

- **`make_bot_og` returns `"rsb_scale"` key** but docstring lists it as `"rsb_scale: The scaling factor"` with no explanation of physical meaning or typical range. RST says "Scaling factor for Pb/Rsb harmonization" — still opaque.

- **`make_pvtw_table` returns `cw_ref` as a 2-element list** (`[cw_usat, cw_sat]`) but the dict key is singular `"cw_ref"` without plural hint. RST/docstring say "Compressibility at reference pressure (1/psi)" — doesn't mention it's a 2-list. Users do `cw_ref * something` and get a list-repeat.

## Low

- **Task hint in prompt mentioned "Fetkovich / Carter-Tracy aquifer"** — no such function exists in simtools. Only `influence_tables` (Van Everdingen-Hurst). Not an API issue, just a spec mismatch.

- **`__all__` lists 15 entries but docstring "Functions" header lists 14** (missing `fit_rel_perm_best` in the module header block, present in `__all__`).

- **`zip_check_sim_deck` returns `None`** when `console_summary=True` and returns `missing` list when False. Two return shapes based on a print-control flag is a code smell — always return `missing` and let the caller ignore.

- **Tier 2 rule cites `deck_check`** (line 10 in `.claude/rules/simtools-module.md`) but the actual function name is `zip_check_sim_deck`. Rule file is stale.

- **`rel_perm_table` RST list (line 448) says `SGFN: Gas / Water table`** — the actual keyword is SGWFN (and SGFN is a different Eclipse keyword for 2-phase gas-water). Documentation naming collision.

- **`ensure_numpy_array` is exported implicitly via `import *`** — not in `__all__` but not prefixed with `_`. Should be `_ensure_numpy_array`.

- **`is_let_physical` returns False silently** when `L/E/T` are such that np.where's denom-guard kicks in; no way for user to distinguish "genuinely unphysical" from "numerically degenerate".

- **`make_vfpprod` silently inserts `0` into `wfr_values` and `gfr_values`** (line 1578–1582). User passes `[0.1, 0.2]` and gets back `[0, 0.1, 0.2]` in the result dict. Mentioned nowhere in docstring.

- **`pref` key in `make_pvtw_table` returns inputs**; `pref` would more conventionally denote a reference constant, not `pi`. Redundant with user-supplied `pi` anyway.

- **Example RST for `rr_solver` returns `0.944, 0.056`** without showing behavior under near-critical `ki~1` (which triggers the correctness bug). Docs promise the method is robust; users assuming this discover otherwise.

---

## Summary (≤150 words)

Counts: 2 Critical, 6 High, 9 Medium, 9 Low.

Top 3:

1. **Inconsistent `export`/return contract across table generators.** VFP returns `"eclipse_string"` in the dict; BOT/PVTW/SWOF/AQUTAB only write to disk — and all file-writing paths write to CWD with no directory control, no `written_path` in the result. Agents can't use these uniformly.

2. **Interactive `input()` prompts in `zip_check_sim_deck` and `ix_extract_problem_cells`** deadlock non-TTY agent use. No non-interactive flag. Also `make_vfp*` silently pads failed cells with 1e-6 BHP and prints (doesn't warn) — decks look valid but simulate wrong.

3. **Parameter/string-enum inconsistencies**: `krtable='SWOF'` vs `well_type='gas'` vs `flo_type='WAT'`; `SGFN`/`SGWFN` naming collision in docs; 7 cryptic method Enums in `make_bot_og` with no inline option list. Discoverability is poor for agents.
