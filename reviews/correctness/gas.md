# Correctness — gas module

## Critical (bug, user gets wrong answer)

- **[gas.py:1945-1946] `_required_concentration` premature exit on small derivative** — When Newton-Raphson derivative `fp` falls below 1e-15, function breaks without computing the next iteration step. This returns an unconverged (incorrect) inhibitor wt% for hydrate suppression. The `w` value remains at its previous iteration rather than being updated via the Newton step (line 1947), violating the convergence contract. Impact: hydrate inhibitor dose may be wrong by orders of magnitude if the cubic tangent slope happens to be very shallow near the solution. Recommended fix: Either (a) compute `w_new` before checking `fp < 1e-15`, then break after assignment, or (b) let the iteration continue naturally and skip the early exit (rely on the convergence check at line 1950).

## High (edge case, silent failure)

None found in this pass.

## Medium (validity range, warning-level)

None found in this pass (correlation validity range warnings are already in place at lines 760-776 for DAK/HY/WYW).

## Low (nitpick)

- **[gas.py:1221] Docstring typo in `gas_bg`** — Parameter docstring says `h2: Molar fraction of Nitrogen` should say `Hydrogen`. No logic error, but misleads users.

