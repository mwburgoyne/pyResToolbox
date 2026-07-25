# API/UX — layer

Scope: `/home/mark/projects/pyResToolbox/pyrestoolbox/layer/layer.py` and `/home/mark/projects/pyResToolbox/pyrestoolbox/docs/layer.rst`.

The public surface is small (5 functions) and mostly clean. Signatures, defaults, and re-exports are reasonable. The main issues are a handful of docstring defects, a signature/doc mismatch, and a couple of magic-negative-sentinel parameters whose UX could be clearer.

## Critical
None.

## High

### 1. `lorenz_2_flow_frac` docstring is internally inconsistent and documents `lorenz` twice with contradictory meaning
`layer.py:181-197`. Two separate `lorenz:` entries appear in the docstring body:

- "lorenz: (0-1) Lorenz hetrogeneity factor"
- "lorenz: Lorenz coefficient (0-1). If B is provided, will ignore this parameter to be more efficient..."

The second entry is the correct precedence logic; the first is an incomplete duplicate. Users and LLM agents reading the docstring will see two conflicting definitions. Also the RST (`layer.rst:176`) describes only the "B overrides lorenz" logic — so RST is cleaner than the code docstring.

### 2. Signature/doc mismatch on `lorenz_2_layers` default for `phi_h_fracs`
`layer.py:225` uses `phi_h_fracs: list = None` (the correct pattern — avoids mutable default), but `docs/layer.rst:212` advertises the signature as `phi_h_fracs = []`. RST should be updated to match, else agents reading the doc will believe passing `None` is invalid.

### 3. `B = -1` sentinel is undiscoverable; type hint claims `float` but semantics are tri-state
`layer.py:182, 224`. Parameter `B: float = -1` is actually a tri-state flag: "<0 means compute from lorenz; >0 means use directly". Type hint `float` gives no signal that negative values are sentinel. Preferred API would be `B: Optional[float] = None` with an explicit check. Low-priority fix but deserves elevation because both `lorenz_2_flow_frac` and `lorenz_2_layers` share the pattern.

## Medium

### 1. `lorenz_2_flow_frac` error message "Must define either B or lorenz parameters" is misleading
`layer.py:201`. The guard triggers only when both `B < 0` and `lorenz < 0`. But valid `lorenz == 0` means homogeneous and should be acceptable — the check `lorenz < 0` rejects only negative values, which is fine, but the error text implies "you forgot one". Should read "`lorenz` must be in [0, 1] unless `B` is provided" or similar, and `lorenz=0` should route to a homogeneous-response short-circuit rather than falling into the solver.

### 2. `lorenz_from_flow_fraction` error messages lack actual values
`layer.py:143-148`. Raises like `"kh_frac must be greater than phih_frac"` without printing `kh_frac` / `phih_frac`. Every other module in the codebase includes the offending values (e.g., `f"Lorenz coefficient must be between 0 and 1, got {lorenz}"` on `layer.py:60` — same module). Inconsistent within a single file.

### 3. Typos/grammar in docstrings and RST
Multiple places:
- "consistant" -> "consistent" (docstring `layer.py:227`, RST `layer.rst:214`)
- "igonore" -> "ignore" (`layer.py:232`)
- "hetrogeneity" -> "heterogeneity" (`layer.py:184`, RST `layer.rst:5`, `layer.rst:7` "Heteogeneity")
- "explictly" -> "explicitly" (`layer.py:229`)
- "phi_h_fracs" in RST example at `layer.rst:149` — the line `"60% of the observed flow comes from 15% of the net thickness"` sits immediately before `.. code-block::` without a blank line, which breaks RST rendering of the directive.

### 4. `lorenz_2_layers`: `nlayers` default of 1 is a footgun
`layer.py:221`. With default `nlayers=1` and no `phi_h_fracs`, the function silently short-circuits to `np.array([k_avg])` — ignores both `lorenz` and `lrnz_method`. No warning. The docstring mentions ">1 needed unless a list of phi_h_fracs is supplied" but there is no error if the user just forgets `nlayers`. Consider raising, or changing default to something non-silent (e.g., `nlayers=10`, or require it).

### 5. `lorenz_2_layers` docstring has `kavg` but signature has `k_avg`
`layer.py:234`. The docstring parameter block writes "kavg:" — a stale name. Signature and RST use `k_avg`. Minor, but trips both grep and LLM consumers.

## Low

### 1. Return type of `lorenz2b` for clamped-boundary cases
`layer.py:62-74`. When `lorenz < 0.000333` or `lorenz > 0.997179...` the function returns hard-coded B limits without any signal. Fine behavior, but no warning is emitted and no docstring note — users looking at `B = 709` won't know it was clamped.

### 2. LinkedIn URL as primary citation
Every function's docstring and the RST quote `https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/` as the canonical reference. LinkedIn URLs are fragile (login walls, link rot). Given this is a published library, a persistent mirror (author's blog, SPE paper if one exists) would be more durable for agents reading docs at runtime.

### 3. `__all__` and re-exports are clean
`layer/__init__.py` is a single-line `from .layer import *`, and `layer.py:34-37` declares `__all__` explicitly. Fine.

### 4. Module docstring in `layer.py` lists 5 functions and matches `__all__` exactly
No discrepancy. Good.

### 5. `bisect_solve` not exposed — correct
Internal dependency; not re-exported. Appropriate.
