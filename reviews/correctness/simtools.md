# Correctness — simtools module

## Critical (bug, user gets wrong answer)

- **simtools.py:1074** — RR solver: `phi_min = 1/(1 - max(ki_hat))` crashes when `max(ki_hat) ≥ 1.0`. Triggered near-critical with K-values ≥ 1. Cascades to line 1078 where `b = 1/(0.25 - phi_min)` also fails if `phi_min ≥ 0.25`.
  **Fix:** Guard `max(ki_hat) < 1.0` (clamp to 0.9999) and add explicit constraint check.

- **simtools.py:1078** — `b` initialization violates bisection/Newton assumptions when `phi_min ≥ 0.25` (negative or undefined `b`). Test case `ki = [1.5, 0.9, 0.1]` triggers immediately.
  **Fix:** Enforce `phi_min < 0.25` or restructure initialization.

- **simtools.py:1679** — VFP metric GFR export scaling: oil-well metric GFR is converted in/out with an extra `*1000` factor, producing a ~1,000,000x magnitude error on export. Line 1544 converts sm3/sm3 → scf/stb, 1630 multiplies by 1000 → Mscf/stb, 1679 multiplies by 1000 again AND converts back to metric.
  **Fix:** Remove the 1000 factor on line 1679 or restructure the chain to be explicit.

## High (edge case, silent failure)

- **simtools.py:1098–1099** — RR solver silently breaks at `max_iter` with no warning; returns unconverged compositions. Regression vs `_lib_vle_engine.rr_solver` which does warn (fixed 2026-04-12).
  **Fix:** Add `warnings.warn("RR solver did not converge in {N_it} iterations")`.

- **simtools.py:350–352** — Saturation normalization returns all zeros when `s_min == s_max` (overlapping endpoints). Produces `kr=0` everywhere even at saturation 1.0 in actual space. Triggered by invalid endpoint definitions (`swcr > 1-sorw`).
  **Fix:** Raise `ValueError`; do not silently collapse.

## Medium (validity / robustness)

- **simtools.py:763–764** — Van Everdingen-Hurst influence-table Python fallback: Bessel-product denominator in `s^1.5 * bessel_product` can approach zero, producing NaN/Inf.
  **Fix:** Guard denominator magnitude; fallback or raise.

- **simtools.py:1680–1684** — VFP metric zero anchoring happens AFTER metric conversion. If user didn't supply 0 in original, anchoring acts on already-converted list → duplicate/malformed axis.
  **Fix:** Anchor zeros BEFORE metric conversion.

## Low (nitpick)

- **RR solver** — References "Nielsen & Lia (2022)" without full citation / equation-number mapping to paper.
- **simtools.py:317** — Jerauld `b = max(b, 1e-10)` silently clamps instead of rejecting `b ≤ 0`.
- **make_vfpprod** — No validation of `alq_values` (accepts negatives / absurd magnitudes).

## Confirmed Correct
Corey, LET, Jerauld formulas match literature; PVTW compressibility unpack handles `[cw_usat, cw_sat]` tuple; `make_bot_og` delegation routes correctly to oil helpers; Bessel formula matches Van Everdingen-Hurst standard.

## Top 3 priorities

1. RR solver division-by-zero when `max(ki_hat) ≥ 1.0` (can crash; also feeds brine VLE).
2. RR solver silent non-convergence (wrong compositions propagate into CO2-solubility calcs).
3. VFP metric GFR oil-well export factor-of-1,000,000 error (simulation deck will be wrong).
