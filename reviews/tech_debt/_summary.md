# Tech Debt Pass — Consolidated Summary

Review date: 2026-04-17. Eleven modules, against prior Tier 3 pass baseline (2026-04-12).

## Severity totals

| Module | High | Med | Low |
|---|---|---|---|
| gas | 3 | 13 | 10 |
| oil | 4 | 8 | 9 |
| brine + plyasunov | 5 | 10 | 9 |
| nodal | 3 | 7 | 8 |
| dca | 2 | 5 | 9 |
| matbal | 2 | 5 | 5 |
| simtools | 5 | 8 | 12 |
| layer | 2 | 8 | 6 |
| plyasunov (standalone) | 3 | 8 | 8 |
| recommend | 3 | 5 | 7 |
| sensitivity | 0 | 5 | 7 |
| **Total** | **32** | **82** | **90** |

## Cross-cutting themes

### 1. Rust parity was under-tested
- **nodal**: Rust has 8 copy-pasted per-method VLP segment loops (~800 lines) — the `_segment_march` scaffold was never propagated to Rust. Every future algorithmic change must be done twice.
- **nodal**: thin Rust-parity tests (8 tests, loose tolerances) would mutually validate the BB `**4` bug as "equal."
- **nodal Rust**: magic numbers never got the Python named-constants treatment (0.0765, 62.4, 5.615, 86400, 141.5, etc. still inlined).
- **brine**: Python + Rust flash run in parallel with no parity test suite — the Sechenov `ks` divergence flagged in correctness is a symptom.
- **matbal**: `objective()` is an ad-hoc Rust try/except (the pattern `@rust_accelerated` was supposed to replace). Lone holdout after Tier 3 nodal migration.

**Recommendation:** add Rust-vs-Python parity test suite with tight tolerances across a representative T/P/composition grid per accelerated module.

### 2. Monoliths survived Tier 3
- **gas/gas.py** — 2195 lines, one file. Hydrate prediction (~500 lines) is a clean carve-out.
- **simtools/simtools.py** — 2064 lines, one file. Same justification as the oil split.
- **brine/_lib_vle_engine.py** — 2726 lines. Split along component/BIP/alpha/flash/API lines.
- **nodal/nodal.py** — 2132 lines. Candidate for a segment-march + IPR + operating-point split.

### 3. Residual duplication
- **simtools**: `rr_solver` duplicated verbatim between `simtools.py:1017` and `brine/_lib_vle_engine.py:1011`. Brine's copy has K=1 guard + convergence warning; simtools' does not (this is *why* correctness flagged two bugs). Delete simtools copy, re-export.
- **simtools**: `ensure_numpy_array` duplicates `shared_fns.convert_to_numpy`; 3 near-identical SWOF/SGOF/SGWFN blocks in `rel_perm_table`; `_format_vfpinj`/`_format_vfpprod` share prologues; 5 independent `open("*.INC","w")` writers with no shared writer.
- **dca**: 6 `_fit_*` helpers share identical R²/residuals boilerplate; two hyperbolic grid-searches structurally identical; `fit_decline` + `fit_decline_cum` duplicate scaffold (~90 lines).
- **gas**: `_h2_method_override + validate_methods` block copy-pasted 10× with asymmetric BNS cross-coupling only in 2 of 10 sites. Single `_resolve_methods()` fixes both duplication and inconsistency.
- **matbal**: Efw/denom/N-mean math triple-maintained (vectorised NumPy + Python scalar fallback + Rust).
- **brine**: Spivey density chain re-implemented inside `brine_props_co2`; `SoreideWhitson._calc_properties` calls `brine_props` 4× for same (p,T).
- **layer**: EXP/LANG dispatch block duplicated 5× across `lorenz2b`, `lorenzfromb`, `lorenz_from_flow_fraction`, `lorenz_2_flow_frac`, `lorenz_2_layers`.
- **plyasunov**: gas-alias maintenance duplicated across 3 dicts (B12 elif, P_IN_MAP, MW_GAS).

### 4. Dead code & unused imports
- **gas**: 22 dead names on import lines (pandas, 9 unused Enums, 12 unused constants); `mws` dead at 864; redundant `zee_arr, _ = convert_to_numpy(zee)` at 1115.
- **oil**: AST-verified dead imports across 5 of 10 sub-files (numpy, psc, oil_sg, oil_pbub, oil_viso, oil_deno, SCF_PER_STB_TO_SM3_PER_SM3, INVBAR_TO_INVPSI); `_cofb_mccain` + `_perrine_co_sat` re-exported but nothing imports them; `oil_bo_mccain_rust` Rust dispatch registered but never called.
- **nodal**: `_gas_density`, `vlp_method`, `class_dic`, 9 constants dead-imported; state-dict has ~11 unused keys; duplicate `_LBFT3_TO_KGM3`.
- **matbal**: dead `List`, `PSI_TO_BAR`.
- **simtools**: 6 unused imports (os, gas, vlp_method, uo_method, 4 unit constants); dead `temp_unit` in both VFP formatters.
- **recommend**: `sg` and `well_type` threaded through two functions, defaulted, documented, never read.
- **plyasunov**: `rho_and_kappa` is dead code that duplicates the work it was meant to avoid.

### 5. Tier-3 work not fully propagated
- **gas, nodal, simtools**: files never split (Tier 3 T3-3 only covered oil).
- **matbal**: never migrated to `@rust_accelerated` (Tier 3 T3-5 covered 8 nodal sites but skipped matbal).
- **gas, nodal, matbal**: still inline some constants that `_constants.py` already defines (141.5/131.5 API conversion, SCF conversions).
- **oil**: `sg_st_gas` and `oil_rs_st` still have inline Valko-McCain matrices the named-constants pass missed.
- **matbal**: missing `validate_pe_inputs` at both public entry points — the only remaining module.

### 6. Coverage gaps
- **matbal**: no metric-units pvt_table test (would have caught the correctness bug).
- **simtools**: `zip_check_sim_deck`, `ix_extract_problem_cells`, `influence_tables`, `make_bot_og`, `make_pvtw_table`, VFP metric mode — all untested.
- **dca**: no Rust-vs-Python parity test for hyperbolic fitters.
- **plyasunov**: no dedicated `test_plyasunov.py` — only exercised indirectly via brine.
- **oil**: `sg_st_gas`, `sgg_wt_avg`, `oil_ja_sg`, `_perrine_co_sat`, `_tables.py` helpers — missing or transitive-only.

### 7. Type hints
Near-universally absent on public signatures across dca, matbal, nodal, oil, simtools, plyasunov. Inconsistent with growing Python project convention.

### 8. Magic numbers that escaped Tier 3
- **dca**: 17 inline magics (b-grid bounds, `maxfev=5000`, `1e-10`, `1e-30`, Duong `0.001`/`500`/`10`).
- **simtools**: `1e-30`, `1e-10`, `1e-6` unnamed.
- **matbal**: `1e-30` ×6, `1e10` ×2, `1e-6` ×2.
- **gas**: `sg_hc` calc has 3 spellings; line 500 uses hardcoded MW literals that don't quite match the named constants.

## Top 10 tech debt items by impact

1. **Add Python↔Rust parity test suites** for nodal, brine, matbal, dca (tight tolerances, representative grid). Would surface BB-exponent bug, Sechenov ks divergence, and future drift automatically.
2. **Propagate `_segment_march` scaffold to Rust nodal** (and propagate named constants).
3. **Split monoliths**: gas.py, simtools.py, _lib_vle_engine.py, nodal.py (follow the oil split template).
4. **Unify the two `rr_solver` copies** (simtools + brine) into one authoritative implementation.
5. **Migrate matbal.objective to `@rust_accelerated`** (last holdout) + add `validate_pe_inputs` to its entry points.
6. **DCA refactor**: shared fit scaffold (R²/residuals/windowing) + named constants + type hints.
7. **Gas `_resolve_methods()` helper**: fixes 10× copy-paste and the asymmetric BNS coupling.
8. **Dead-code cleanup**: oil (`_cofb_mccain` re-export, `oil_bo_mccain_rust` dispatch, 5-file import leak); gas (22 unused names); nodal (state-dict + 9 constants); recommend (`sg`/`well_type` dead params).
9. **Layer EXP/LANG dispatch** → 2 private helpers reduces 5 copies to 1.
10. **Plyasunov gas-alias consolidation** + add dedicated test file.

## Sensitivity is the cleanest
Zero High findings; only polish-level items.
