# Correctness — nodal module

## Critical (bug, user gets wrong answer)

- **nodal.py:1543 AND src/vlp/holdup_bb.rs:111** — Beggs-Brill two-phase friction S-factor polynomial uses `ln_y ** 4` instead of `ln_y ** 3`.
  The polynomial should be `c0 + c1*ln(y) + c2*ln(y)^2 + c3*ln(y)^3` (cubic). Current code uses a quartic term, causing divergence at high liquid holdup (high y).
  **Impact:** Incorrect friction correction for ALL Beggs-Brill gas/oil calculations across all inclinations. Affects BHP predictions.
  **Fix:** Change both instances from `** 4` to `** 3`. Update both Python and Rust together (parity requirement).

## High (edge case, silent failure)

- **nodal.py:2062** — Operating-point IPR interpolation uses `np.interp()` on reversed IPR curve without bounds checking. If the bisection proposes a rate outside IPR bounds, extrapolation can produce unphysical Pwf. Mitigated by clamps at lines 2067–2068 but fragile.
  **Fix:** Add explicit clamping before the interp call and raise/warn when clipping occurs.

- **nodal.py:1163–1165** — Gray holdup logarithm: argument `(a + b*n2)/n1` can approach 0 at high Froude numbers, causing `ln(~0) → -inf` before downstream clamping kicks in.
  **Fix:** Floor the argument (`if argument < 1e-10: f_e = 0.0`) before the log.

## Medium (validity / robustness)

- **nodal.py:1569–1579** — BB TRANSITION flow pattern recomputes horizontal holdup twice; functionally correct but duplicative.
- **Condensate dropout linear model** — known Tier 4 limitation, not re-flagged.
- **nodal.py:975–976** — Serghides Fanning factor clamps `Re ≥ 2100`, preventing true laminar friction. Slightly unconservative; document as intentional.

## Low (nitpick)

- **nodal.py:767–782** — IFT correlations have no API/temperature range validation (mitigated by floor at 1.0 dyne/cm). Consider optional warnings.
- **nodal.py:1067–1077** — WG Chisholm C boundary is discontinuous at `Re = 2100` — matches industry standard, no action.

## Top 3 priorities

1. **Critical BB exponent bug** — fix both Python and Rust, add regression test against a literature benchmark.
2. **IPR extrapolation** — tighten bounds handling in operating-point solver.
3. **Gray log instability** — add floor guard.
