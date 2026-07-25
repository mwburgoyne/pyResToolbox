# Correctness Pass — Consolidated Summary

Review date: 2026-04-17. Eleven modules reviewed in parallel against the 2026-04-12 cold-eyes baseline.

## Severity totals

| Module | Crit | High | Med | Low |
|---|---|---|---|---|
| gas | 1 | 0 | 0 | 1 |
| oil | 3 | 2 | 2 | 3 |
| brine + plyasunov (combined) | 3 | 3 | 7 | 3 |
| nodal | 1 | 2 | 3 | 2 |
| dca | 3 | 2 | 4 | 0 |
| matbal | 1 | 4 | 2 | 3 |
| simtools | 3 | 2 | 2 | 3 |
| plyasunov (standalone) | 1 | 2 | 1 | 2 |
| layer | 0 | 0 | 0 | 0 |
| recommend | 0 | 0 | 1 | 2 |
| sensitivity | 0 | 2 | 3 | 2 |
| **Total** | **16** | **19** | **25** | **21** |

Note: plyasunov counted both inside brine (as its call-site impact) and standalone (for IAPWS-IF97 specifics).

## Top 10 critical issues (fix first)

1. **[nodal.py:1543 + src/vlp/holdup_bb.rs:111]** Beggs-Brill friction S-factor polynomial uses `ln_y**4` instead of `ln_y**3`. Affects every BB gas/oil BHP calculation. Fix both Python and Rust together.
2. **[brine/brine.py:1267-1273 garciaDensity]** Garcia density singularity at `xCO2 → 1` is unguarded — regression of a prior fix. Silent NaN/zero density in CO2 storage scenarios.
3. **[brine/_lib_vle_engine.py Python vs src/vle/flash.rs]** Python/Rust Sechenov ks divergence for CO2/H2S. Rust silently uses S&W Eq 8 fallback, so `framework='proposed'` with Rust accelerator is wrong by 5–10% in K-values.
4. **[oil/_correlations.py:336]** Rs_velarde NaN when pb ≤ psc (atmospheric bubble points — `0/0`).
5. **[oil/_separator.py:30-40]** `sg_evolved_gas` divides by `p` and `sg_sp` without guards → silent NaN on `p=0` or `sg_sp=0`.
6. **[oil/_density.py:116]** `rhoa` polynomial can zero-out during iteration → silent NaN division.
7. **[simtools.py:1074, 1078]** RR solver `phi_min = 1/(1 - max(ki_hat))` crashes when `max(ki_hat) ≥ 1.0`; cascades to `b` initialization.
8. **[simtools.py:1679]** Oil-well VFP metric GFR export has factor-of-~1,000,000 magnitude error (double `*1000` + double unit conversion).
9. **[matbal P/Z metric pvt_table]** Gas P/Z with `metric=True` and tabulated PVT converts pressures to psia but still uses original bars to interpolate (lines 250, 256). No metric pvt_table test exists.
10. **[plyasunov/iapws_if97.py:137]** Isothermal compressibility formula appears wrong — returns `-gpp/(gp*P_STAR)` but per IAPWS standard should be `-gpp/(gp**2*P_STAR)`. Cascades into V2_inf and SoreideWhitson Garcia mixing.

## Notable High-severity clusters

- **Rust/Python divergence** is a recurring theme (BB exponent, Sechenov ks, Rust fallback swallowed silently in matbal). Suggests automated parity tests are missing.
- **Silent non-convergence** — multiple solvers break at max_iter without warning (simtools RR, CO2_Brine_Mixture Spycher). Was fixed for `_lib_vle_engine.rr_solver` in the 2026-04-12 pass; not propagated to siblings.
- **NaN/Inf propagation in sensitivity tornado** — base-case NaN breaks sort order silently.
- **DCA edge cases** — `dt=0` ZeroDivisionError, negative uptime accepted, non-monotone Np accepted.

## Clean modules
- **layer** — no findings across 12 tests. Mathematical verification clean.
- **recommend** — decision-tree logic sound, method strings valid.

## Next steps
Proceed to Pass 2 (API/UX) and Pass 3 (Tech Debt), then present a consolidated prioritized fix list.
