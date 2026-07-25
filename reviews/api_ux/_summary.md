# API/UX Pass — Consolidated Summary

Review date: 2026-04-17. Eleven modules reviewed in parallel.

## Severity totals

| Module | Crit | High | Med | Low |
|---|---|---|---|---|
| gas | 0 | 6 | 16 | 9 |
| oil | 2 | 6 | 7 | 9 |
| brine | 1 | 5 | 7 | 9 |
| nodal | 2 | 8 | 9 | 8 |
| dca | 2 | 5 | 7 | 7 |
| matbal | 2 | 5 | 6 | 7 |
| simtools | 2 | 6 | 9 | 9 |
| layer | 0 | 3 | 5 | 5 |
| plyasunov | 0 | 2 | 6 | 5 |
| recommend | 0 | 2 | 6 | 6 |
| sensitivity | 0 | 2 | 8 | 6 |
| **Total** | **11** | **50** | **86** | **80** |

## Cross-cutting themes

### 1. Docs/code drift
- **gas**: BUR vs BNS name confusion between docstrings (BNS) and RST (BUR); `gas_grad2sg` SG bounds differ across docstring (0.55–3.0), RST (0.55–1.75), code (0.07–3.0).
- **oil**: RST uses `sg_sto` where code expects `sg_o` — copy-paste from docs fails.
- **oil**: BOT dict key `rsb_scale` vs everywhere-else `rsb_frac`; `vis_frac` not in RST dict-return table.
- **simtools**: docstring/RST say `SGFN` but code expects `SGWFN`.
- **simtools**: Tier 2 rule still references old `deck_check` name.

### 2. Invalid-method errors don't list valid options
- **nodal**: `"An incorrect method was specified: 'INVALID' for 'vlpmethod'"` — no list of HB/WG/GRAY/BB.
- **simtools**: 7 cryptic BOT method enums (EXPLT, PMC, VELAR, SWMH, MCAIN…) with no inline option list.
- Pattern appears module-wide — agent-facing UX is bad here.

### 3. Silent/unvalidated parameters
- **nodal**: `well_type` is never validated; typos silently take the oil branch.
- **gas**: `gas_ug` silently discards `zee` on length mismatch.
- **gas**: `GasPVT` silently drops user-supplied `tc`/`pc`.
- **nodal**: `gas_pvt` param accepted but "reserved for future use" — dead, masks effects.
- **brine**: `ppm` vs `wt` alias resolution silently picks one when both passed.
- **recommend**: `sg`, `well_type` parameters documented but unused by decision logic.

### 4. Hidden configuration knobs
- **brine** (Critical): SoreideWhitson hardcodes `framework='proposed'` and `salinity_method='gamma_phi'` — no constructor override. RST doesn't mention the knobs. Users can't select alternative frameworks even though the engine supports them.
- **nodal**: `injection=True` only honored by direct `fbhp`; `outflow_curve` / `operating_point` silently ignore it.

### 5. Output shape inconsistency
- **brine**: `CO2_Brine_Mixture.Rs` is float; `SoreideWhitson.Rs` is dict → non-portable code.
- **brine**: `SoreideWhitson.y` documented as "gas phase" but is actually input-normalized dry-gas; `CO2_Brine_Mixture.y` is true flash vapor (same name, different meanings).
- **nodal**: `outflow_curve` → `'rates'`; `ipr_curve` → `'rate'`; `operating_point` → `'rate'` (inconsistent keys).
- **dca**: `DeclineResult` uses sentinel `0.0` for unused fields (Arps vs Duong); `RatioResult.a/b/c` positionally overloaded with method-dependent meanings.
- **simtools**: Table generators return/write inconsistently — `make_vfpprod` returns `eclipse_string` in dict; `make_bot_og`/`make_pvtw_table`/`rel_perm_table` only side-effect write `.INC` to CWD with no `written_path` in return.

### 6. Missing integrations
- **matbal** (Critical): no aquifer plug-in. `oil_matbal` has no `We` param; `gas_matbal` requires precomputed `We`. `simtools.influence_tables` exists but isn't wired; no Fetkovich/Carter-Tracy helper.
- **matbal**: no plotting hooks for Cole plot, P/Z vs Gp, Havlena-Odeh, drive-index chart.
- **dca**: no DataFrame convenience I/O despite DCA being a DataFrame-heavy domain.

### 7. Naming/default drift
- **oil**: `pbmethod` default VALMC in some fns, VELAR in harmonize/BOT/`from_harmonize` — chains silently use two Pbs.
- **simtools**: parameter casing mixed (`krtable='SWOF'` upper vs `well_type='gas'` lower vs `flo_type='WAT'` upper).
- **nodal**: `api=45` default fine for condensate but wrong for oil — silent bad answers.
- **brine**: `metric` default differs across the three models (one False, two True) — breaks Tier 1 convention.

### 8. Non-interactive/agent-hostile UX
- **simtools**: `zip_check_sim_deck` and `ix_extract_problem_cells` prompt via `input()` — dead-locks non-TTY callers. No `non_interactive` flag.
- **dca**: no `__repr__` guard on array fields; printing a result floods the terminal.
- **nodal**: error messages echo post-conversion values (`-0.3048 ft` when user typed `-1 m`).

## Top 10 fixes (highest user impact)

1. **brine** — expose `framework` and `salinity_method` on `SoreideWhitson.__init__`; document in RST.
2. **nodal** — validate `vlp_method` and `well_type` inputs; list valid options in the error message.
3. **matbal** — wire up aquifer helpers (Fetkovich + Carter-Tracy + VEH) into `gas_matbal` and `oil_matbal`.
4. **oil** — fix RST `sg_sto` → `sg_o`; align BOT `rsb_scale` → `rsb_frac`; add `vis_frac` to RST dict table.
5. **gas** — canonicalise BUR/BNS naming across code + RST; align `gas_grad2sg` SG bounds in all three places.
6. **simtools** — stop using `input()` in non-interactive paths; uniform return contract for table generators with `written_path` + `eclipse_string` + explicit `out_dir`.
7. **brine** — unify `Rs` return type (dict); rename `SoreideWhitson.y` to avoid meaning clash.
8. **nodal** — return pressure-traverse output + implement `fthp` (THP-from-BHP); expose `converged` flag from `operating_point`.
9. **dca** — add DataFrame input/output convenience; fix `forecast()` t=0 omission; add `__repr__` guard.
10. **matbal** — add plotting hooks (Cole, P/Z, HO straight-line, drive indices) or return figs for user.

## Clean-ish modules
- **layer**, **plyasunov**, **recommend**, **sensitivity** — no Critical findings. Mostly polish-level items.
