# pyResToolbox "ultrareview" — 2026-04-17

Three-pass automated review across 11 modules. Agents ran in parallel per pass.

## Structure

```
reviews/
├── correctness/   # Pass 1: numerical bugs, edge cases, validity ranges
│   ├── _summary.md
│   └── <module>.md × 11
├── api_ux/        # Pass 2: signatures, defaults, docs/code drift, error msgs
│   ├── _summary.md
│   └── <module>.md × 11
└── tech_debt/     # Pass 3: duplication, dead code, split candidates, coverage gaps
    ├── _summary.md
    └── <module>.md × 11
```

## Severity totals (all passes combined)

|  | Critical | High | Medium | Low | Total |
|---|---|---|---|---|---|
| Correctness | 16 | 19 | 25 | 21 | 81 |
| API/UX | 11 | 50 | 86 | 80 | 227 |
| Tech Debt | — | 32 | 82 | 90 | 204 |

## The 10 most important fixes

Ranked across all three passes. "Impact" is a subjective reviewer call, not a severity level.

1. **Beggs-Brill friction S-factor exponent bug** — `nodal.py:1543` + `src/vlp/holdup_bb.rs:111` use `ln_y**4` instead of `ln_y**3`. Affects every BB gas/oil BHP calculation; both Python and Rust are wrong together. *Correctness / Critical.*
2. **Garcia density singularity regression** — `brine.py:1267-1273` guard for `xCO2 → 1` was present in the 2026-04-12 review but is not in current code. Silent NaN/zero density in CO2 storage. *Correctness / Critical.*
3. **Python/Rust Sechenov ks divergence** — Rust `flash_tp` silently uses S&W Eq 8 for CO2/H2S when Python routes to Dubessy/Akinfiev. `framework='proposed'` + Rust accelerator = wrong K-values by 5–10%. *Correctness / Critical.*
4. **Oil NaN bugs** — Rs_velarde `0/0` at atmospheric pb (`oil/_correlations.py:336`), `sg_evolved_gas` divide-by-zero (`_separator.py:30-40`), `rhoa` polynomial zero-out (`_density.py:116`). All silent. *Correctness / Critical.*
5. **simtools RR division-by-zero + oil-well VFP 1M× error** — `simtools.py:1074, 1078` crash when `max(ki_hat) ≥ 1`; `simtools.py:1679` double-applies `*1000` on metric GFR export. *Correctness / Critical.*
6. **matbal metric pvt_table bug** — `metric=True` converts pressures to psia but still interpolates in bars (lines 250, 256). No metric-pvt_table test. *Correctness / Critical.*
7. **Hidden brine framework knobs** — `SoreideWhitson` hardcodes `framework='proposed'` and `salinity_method='gamma_phi'`; users can't pick the alternatives even though the engine supports them. *API/UX / Critical.*
8. **matbal aquifer integration missing** — `oil_matbal` has no `We` param; `gas_matbal` requires precomputed `We`; VEH exists in simtools but isn't wired. Canonical water-drive workflow unsupported. *API/UX / Critical.*
9. **Nodal missing fthp + pressure-traverse output + unvalidated well_type** — listed but not implemented. *API/UX / Critical.*
10. **Python/Rust parity test gap** — multiple modules (nodal, brine, matbal, dca) have thin-or-missing parity tests. Items 1 and 3 would have been caught automatically with tight-tolerance grids. *Tech Debt / High.*

## Cross-cutting themes

- **Silent failure patterns**: non-convergence swallowed by RR solvers (simtools, CO2_Brine_Mixture), Rust fallback with no warning (matbal), NaN propagation in sensitivity tornado, invalid-method error messages that don't list valid options (nodal, simtools).
- **Docs/code drift**: BUR vs BNS (gas), `sg_sto` vs `sg_o` (oil), `SGFN` vs `SGWFN` (simtools), `rsb_scale` vs `rsb_frac` (oil BOT), 3 different `gas_grad2sg` SG bounds.
- **Tier-3 work not fully propagated**: gas/simtools/brine/_lib_vle_engine monoliths still not split; matbal is the last `@rust_accelerated` holdout; Rust side of nodal didn't get the named-constants treatment.
- **Agent-hostile UX**: interactive `input()` in simtools, invalid-method errors without option lists, non-guarded `__repr__` in dca that floods stdout.

## Coverage notes / caveats

- Agents referenced the 2026-04-12 cold-eyes review to avoid re-flagging already-fixed items, but this isn't a guarantee of zero overlap.
- Two read-only (Explore) agents in Pass 1 couldn't write their files; their full output was captured via the handoff message and written manually (brine, sensitivity, nodal, simtools).
- Pass 1 total counts include `plyasunov` both as standalone and as part of brine (the dual count overstates by ~3–4).

## Proposed next steps

Before implementing: suggest the user scan the `_summary.md` files and Top-10, and pick a scope.

Reasonable slicing of the follow-up work:

- **Release-blocking (numerical bugs)**: Top-10 items 1–6. One focused PR fixing each bug with a regression test.
- **Next minor release (v3.2)**: API/UX Critical + matbal aquifer integration + hidden brine knobs.
- **Next Tier-3 maintenance pass**: monolith splits, Rust parity test suites, dead-code sweep, type-hint backfill.

## Output files

All findings per module are in their respective subfolder. The `_summary.md` in each pass has per-module counts, cross-cutting themes, and top-N lists. This README is the tl;dr.
