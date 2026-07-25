# Consolidated Findings — 2026-04-17 (v2 after v3.1.5 + deeper verification)

Synthesis across Correctness / API-UX / Tech Debt passes, **narrowed after per-claim verification against source and published literature**. Earlier reviewer claims that turned out to be false positives are recorded at the bottom so future passes don't re-flag them.

**Status as of 2026-04-17:** v3.1.5 shipped (WS-1 trimmed + WS-9b). Remaining scope split into v3.2.0 / v3.3.0 / v3.4.0 below.

---

## Part 1 — Release plan

Three releases, each a coherent theme, each reviewable as a handful of focused PRs rather than one megabundle.

### v3.2.0 — Correctness + tiny ergonomics (days)
Low-risk, high-value. Bundles the remaining **verified** correctness bugs with small ergonomic helpers.

| WS | Content | Risk |
|---|---|---|
| **WS-10 bugs** | `duong_cum` negative cum for `t < 0.001` (line 329 linspace inversion); `forecast(dt=0)` ZeroDivisionError; `forecast` missing `uptime ∈ (0,1]` validation | Low — DCA-only, 3 localised fixes |
| **WS-2 dedup** | Delete `simtools.rr_solver` duplicate, re-export from `brine._lib_vle_engine`. Delete `simtools.ensure_numpy_array`, use `shared_fns.convert_to_numpy`. K=1 regression test. | Low — pure dedup |
| **WS-7 Sechenov guard** | Rust `ks` only implements S&W Eq 8; Python `framework='proposed'` uses Dubessy (CO2) / Akinfiev (H2S) specialised `ks`. Impact is **zero today** because Python precomputes γ and passes it to Rust, but add a guard + docstring so the divergence can't silently bite if someone uses Rust `ks` standalone later. | Low — no behaviour change, doc + assertion only |
| **WS-9 converged contract** | Define a simple contract: solvers expose `.converged` attribute / dict key. Apply to `CO2_Brine_Mixture` (add `.converged`). Defer SoreideWhitson and `nodal.operating_point` to v3.3 to keep this PR narrow. | Low — additive |
| **WS-14 tiny** | `sensitivity`: NaN guard on `base_result`, validate `lo <= hi`. `layer`: EXP/LANG dispatch → helpers (5 copies → 1). `recommend`: remove unused `sg`/`well_type` params. Verify each before patching. | Low — isolated |

### v3.3.0 — API completion + non-breaking ergonomics (weeks)
Coherent around "finish the API surface without changing behaviour."

- **WS-5** nodal API: `fthp` (THP-from-BHP), pressure-traverse output from `fbhp`, `vlp_method` validation via the new helper, `converged` flag on `operating_point`, dict-key unification across curve/IPR/op-point, honour `injection=True` in `outflow_curve`/`operating_point`, drop dead `gas_pvt` param, judgment call on `api=45` default, floor guard on Gray holdup log arg.
- **WS-6** oil hygiene: RST sync (`sg_sto` → `sg_o`, `rsb_scale` → `rsb_frac`, add `vis_frac`), unify `pbmethod` default (VALMC vs VELAR inconsistency across `oil_rs`/`oil_co`/`oil_harmonize`/`make_bot_og`), dead imports across 5 sub-files, drop never-called `oil_bo_mccain_rust` dispatch, move `_cofb_mccain` to `_compressibility.py`, extract remaining Valko-McCain coefficient matrices to `_constants.py`.
- **WS-7 API** (non-breaking): expose `framework` and `salinity_method` on `SoreideWhitson.__init__`, document in RST, validate combos (reject/warn `proposed`+`embedded`, warn `sw_original`+`gamma_phi` double-count), error on both `ppm` and `wt`, unify `metric` default across all three brine models.
- **WS-8** gas hygiene: align `gas_grad2sg` SG bounds across docstring/RST/code, introduce `_resolve_methods()` helper (collapses 10× copy-paste + fixes asymmetric BNS coupling in 2 of 10 sites), don't silently discard user `zee` in `gas_ug` — warn + recompute; don't silently drop `tc`/`pc` in `GasPVT` — accept or error; remove dead imports (`pandas as pd`, 9 unused Enums, 12 unused constants, `mws`).
- **WS-10 polish**: judgment call #3 on hyperbolic b-bounds (`[0,1]` vs `[0,2]`), DataFrame I/O convenience, type hints, extract shared fit scaffold (~90 lines duplicated across `_fit_*`), 17 magic numbers → named constants, fix `Qcum`/`Np` vs `t`/`q` PascalCase inconsistency, fix `forecast` t=0 omission (`np.arange(dt, ...)` skips t=0 — design decision to revisit).
- **WS-9 parity harness** (shared Rust-vs-Python test helper). Feeds into v3.4.

### v3.4.0 — Structural + features + breaking (when you have review capacity)
Monolith splits and feature additions that shouldn't ride with behaviour changes.

- **WS-3 matbal overhaul**: wire aquifer helpers (Van Everdingen-Hurst via `simtools.influence_tables`, add Fetkovich + Carter-Tracy) into `gas_matbal`/`oil_matbal`; add plotting hooks (Cole, P/Z vs Gp, Havlena-Odeh, drive-index stacked); add `validate_pe_inputs` at public entry points; migrate `objective()` to `@rust_accelerated`; consolidate Efw/denom/N-mean triplicate math; dead code removal. *Note: the "metric pvt_table psia/bar bug" and "Havlena-Odeh (1+m)*Efw" originally flagged in WS-3 are verified not-a-bug and removed from scope.*
- **WS-4 nodal Rust sync**: propagate `_segment_march` scaffold to Rust (removes ~800 lines Rust duplication); propagate named constants to Rust; dead-code (`_gas_density`, `vlp_method`, `class_dic`, 9 constants, 11 state-dict keys, duplicate `_LBFT3_TO_KGM3`). *The "BB exponent bug" originally flagged in WS-4 is verified not-a-bug (Beggs & Brill 1973 Eq 16 confirms `**4`) — keep the Rust-sync scope, drop the bug claim.*
- **WS-8 hydrate split**: carve ~500 lines of hydrate code out of `gas.py` (2195 lines) along the clean boundary. Keep public API identical.
- **WS-11 simtools structural**: split `simtools.py` (2064 lines) into rel-perm / VFP / BOT / aquifer / RR; unify SWOF/SGOF/SGWFN blocks in `rel_perm_table` (3 near-identical copies); unify `_format_vfpinj`/`_format_vfpprod` prologues; shared `.INC` writer helper with `out_dir` + returned `written_path`; uniform return contract for table generators (`eclipse_string` + `written_path` in dict); fix `SGFN`/`SGWFN` docstring/RST/code mismatch; 6 unused imports, `temp_unit` dead in both VFP formatters.
- **WS-12 brine VLE engine structural**: split `_lib_vle_engine.py` (2726 lines) along component / BIP / alpha / flash / API boundaries; extract Spivey density-chain helper (re-implemented in `brine_props_co2`); cache `brine_props` calls in `SoreideWhitson._calc_properties` (currently called 4× for the same `(p,T)`).
- **WS-13 plyasunov cleanup**: consolidate gas aliases across 3 dicts (B12 elif, P_IN_MAP, MW_GAS); add dedicated `test_plyasunov.py`; document Region 1 validity (273–623 K); dead `rho_and_kappa`. *The "compressibility formula bug" claim is verified not-a-bug — formula matches the official IAPWS `iapws` library.*
- **Breaking (if user approves)**: `SoreideWhitson.Rs` return type change (float → dict) to unify with `CO2_Brine_Mixture`; `SoreideWhitson.y` rename. Deprecation shim for one minor cycle, then remove. See judgment calls #1 and #2.
- **Type hints backfill** across all modules (judgment call #10).

---

## Part 2 — Judgment calls (still open; need user direction)

1. **`SoreideWhitson.Rs` return type** — float → dict would unify with `CO2_Brine_Mixture`. Add deprecation shim (`.Rs_total` attribute) in v3.3, remove in v3.4? Or hold entirely for v4.0 major bump?
2. **`SoreideWhitson.y` rename** — same deprecation question.
3. **DCA hyperbolic b-bounds** — `[0, 1]` (conservative) vs `[0, 2]` (literature permits for ultra-steep declines). Pick one; v3.3 scope.
4. ~~Matbal Havlena-Odeh `(1+m)*Efw` canonical vs agent claim~~ — **deferred**; needs Dake reference check before v3.4.
5. **Python/Rust Sechenov ks long-term** — v3.2 ships the short path (Rust guard + document). Long-term: port Dubessy/Akinfiev to Rust (weeks of work, low payoff since Python already precomputes γ). Skip unless Rust `ks` gets a new caller.
6. **Condensate dropout** — prototype lives at `/home/mark/projects/CondensateDropout`. Integrate or keep linear-model warning? Not blocking any release.
7. **Nodal `api=45` default** — fine for condensate, wrong for oil. Make required, or keep default + warn? v3.3 scope.
8. ~~Plyasunov isothermal compressibility `gp**2` claim~~ — **closed**; verified not-a-bug against official IAPWS library.
9. **Monolith splits (v3.4)** — do them standalone (pure mechanical, easy to review) or bundle with behaviour changes? Recommendation: standalone, which is what v3.4 plan already does.
10. **Type hints backfill** — one big PR in v3.4 or module-by-module per WS? Recommendation: per-WS as each lands.

---

## Part 3 — Verified not-a-bug (agent false positives)

Recording here so future reviews don't re-flag. All verified against source + published references.

**From v3.1.5 implementation (2026-04-17 first pass):**
- **Beggs-Brill friction S-factor `ln_y**4`** (`nodal.py:1543`, `src/vlp/holdup_bb.rs:111`). Correct per Beggs & Brill 1973 Eq 16 and Brill & Mukherjee SPE Monograph 17 Eq 3.76.
- **simtools VFP metric GFR "1M× error"** (`simtools.py:1679`). Conversion chain is self-consistent — line 1544 `/1000.0` was missed by the reviewer; line 1630 mutates a loop-local only.
- **matbal metric pvt_table** (`matbal.py:250, 256`). `np.interp(p, pvt_p, ...)` consistent (both bars when metric=True); `p_psia`/`pvt_p_psia` used only in Z↔Bg EOS at 252/258 (consistent in psi). No unit mixing.
- **oil `_density.py:116` rhoa polynomial "zero-out"**. Only triggers at `sg_sp ≈ 0.1` (lighter than H₂, non-physical). Defensive programming suggestion, not a bug.

**From v3.2.0 pre-planning verification (2026-04-17 second pass):**
- **Brine alpha function `Tr > 1`** (`_lib_vle_engine.py:127-141, 144-160`; `src/vle/alpha.rs:14-32`). Both `alpha_water_soreide` and `alpha_water_mc3` clamp `Tr_safe = max(Tr, 0.01)` and return positive, physically sensible values across the full `Tr ∈ [0.5, 2.0]` range. Rust matches.
- **Brine Rachford-Rice K-value clipping order** (`_lib_vle_engine.py:1041, 1132-1137, 2143, 2153`). The reviewer's "clip before single-phase check" claim is moot — K=1 perturbation at line 1041 inside `rr_solver` shields the solver regardless of where final clipping happens.
- **Gas `_required_concentration` premature-exit** (`gas.py:1942-1955`). Loop updates `w` in both converged and non-converged branches; round-trip test `test_hydrate_required_concentration_round_trip` passes to machine precision (<1e-16) across depression values 0.1–500.
- **Plyasunov isothermal compressibility `-gpp/(gp*P_STAR)` vs `gp**2`** (`plyasunov/iapws_if97.py:137`). Formula is correct per IAPWS-IF97 Revised Release 1997 Eq 7 and matches the official `iapws` Python library. Reviewer misread dimensionless vs dimensional derivatives.

**Still a real divergence but low-impact:**
- **Sechenov Rust ks for `framework='proposed'`** — Rust only has S&W Eq 8; Python 'proposed' uses Dubessy (CO2) + Akinfiev (H2S). Impact is **zero today** because Python precomputes γ and passes to Rust. v3.2 adds a guard + doc to prevent future misuse (no behaviour change).

---

## Part 4 — What we're NOT touching

- **Low-severity nitpicks** in per-module reviews — skip unless adjacent work makes them free.
- **Aesthetic decisions** (DCA `DeclineResult` sentinel `0.0` pattern) — depends on user preference. Only act if you have one.
- **Speculative correctness claims** that haven't been verified. Current open verifications: matbal Havlena-Odeh `(1+m)*Efw`.
- **Condensate dropout "improve"** — scoped project at `/home/mark/projects/CondensateDropout`, not a review finding.

---

## Summary counts

| | v3.1.5 | v3.2.0 | v3.3.0 | v3.4.0 |
|---|---|---|---|---|
| WS items | WS-1 (trimmed) + WS-9b | WS-10 bugs, WS-2, WS-7 guard, WS-9, WS-14 | WS-5, WS-6, WS-7 API, WS-8 hygiene, WS-10 polish | WS-3, WS-4, WS-8 split, WS-11, WS-12, WS-13, breaking |
| Risk | Shipped ✅ | Low | Medium | Medium-high (structural churn) |
| Est PRs | 1 | 1 coherent | 4–5 per-module | 6–7 per-module + breaking |

Verification has cut 4 "critical" claims to zero across v3.2 and v3.3 scope since the original plan. Net remaining correctness bugs in the library: **3** (all DCA, all in v3.2.0).
