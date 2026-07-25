# Correctness — oil module

## Critical (bug, user gets wrong answer)

- **_correlations.py:336** — Rs_velarde division by zero when pb ≤ psc. After clamping both p and pb to max(psc, x) at lines 320-321, if input pb ≤ 14.696, then pb becomes exactly psc, causing `(p-psc)/(pb-psc) = 0/0 = NaN`. Impact: Rs calculation silently returns NaN when bubble point is at or below atmospheric. Recommended fix: Guard after clamping: `if pb <= psc: raise ValueError("pb must be > psc")` or skip Velarde and return error.

- **_separator.py:30-40** — sg_evolved_gas division by zero. Line 29-39 computes `one_on_sgr = a[1]/p + a[2]/p**2 + ... + a[8]/sg_sp + ...` without guards on p or sg_sp. If p=0 or sg_sp=0 (edge inputs that pass validate_pe_inputs), division by zero occurs. Impact: Silent NaN on edge-case inputs. Recommended fix: Add `if p <= psc or sg_sp <= 0: raise ValueError(...)` before computation.

- **_density.py:116** — Potential division by zero in rhoa. Line 107-114 computes `rhoa = a[0] + a[1]*sg_sp + a[2]*sg_sp*rho_po + ... ` where a[0] = -49.8930 (large negative). The polynomial can zero out for certain rho_po values, then line 116 divides by rhoa: `new_rho_po = (...) / (_SWMH_MASS_DENOM + rs*sg_sp/rhoa)`. When rhoa=0, division by zero. Impact: Silent NaN in density iteration when converging on rho_po. Recommended fix: Guard in iteration loop: `if abs(rhoa) < 1e-10: raise ValueError(...)` after computing rhoa.

## High (edge case, silent failure)

- **_correlations.py:222** — Velarde Rs inverse uses complex arithmetic without explicit handling. Line 222 in rsbub_velarde: `rsb = (-10**(-x)*...*(...)**(1/_VEL_PB_POW) - _VEL_PB_OFFSET)**(1/_VEL_RS_EXP)`. The inner term can be negative, raising a negative base to fractional exponent 1/_VEL_RS_EXP = 12.28, producing complex numbers. Line 227 returns rsb directly (not wrapped in get_real_part). Impact: Complex-valued Rs returned to user when correlation extrapolates to low pressures. Recommended fix: Wrap final rsb: `return get_real_part(rsb)` at line 227.

- **_rate.py:94-106** — Darcy productivity index silently fails if uo or bo are zero. Line 94-96: `J = _DARCY_RADIAL * k * h / (uo * bo * ...)`. No check that uo > 0 and bo > 0. If PVT methods return 0 (rare edge case), user gets Inf/NaN without diagnostic. Impact: Silent Inf in productivity index. Recommended fix: Add guard in oil_rate_radial before J calculation: `if uo <= 0 or bo <= 0: raise ValueError("uo and bo must be positive")`.

- **_harmonize.py:103** — Newton iteration without backtracking. Line 103-111 solves `pb_target = oil_pbub(rsb_new)` via fixed-point iteration `rsb_new = (pb/pbcalc) * rsb_old`. No damping (e.g., 0.5*rsb_old + 0.5*rsb_new), so oscillation is possible for certain pb/rsb combinations. Line 112 raises RuntimeError after 100 iterations, but intermediate NaN/Inf can silently propagate if iteration diverges before 100 steps. Impact: Possible divergence without warning. Recommended fix: Add damping factor (e.g., 0.7) in line 104 or check for NaN/Inf in loop.

## Medium (validity range, warning-level)

- **_correlations.py:75-78, 80-83** — Standing and Valko-McCain Pb correlations issue only warnings outside calibration range, no re-validation on output. Standing: API [16.5, 63.8], T [100, 258]°F; Valko-McCain: API [6, 56.8], T [78, 330]°F. If user inputs violate these, warning issued but Pb correlation executes anyway. No guarantee output Pb is physically reasonable. Impact: User may not notice warning and use invalid Pb. Recommended fix: (non-blocking) Add runtime output checks, e.g., `if pbcalc < 0: warnings.warn("negative Pb, likely outside correlation range")`.

- **_harmonize.py:120-126** — Viscosity scaling vis_frac may become very large or very small without bounds. Line 126: `vis_frac = uo_target / uo_corr`. If uo_corr is very small (e.g., 0.01 cP floor), vis_frac could be huge. No clamping to physically reasonable range (e.g., [0.5, 2.0]). Impact: Extreme viscosity scaling can distort black-oil model. Recommended fix: (optional) Warn if vis_frac < 0.5 or > 2.0.

## Low (nitpick)

- **_correlations.py:106-108** — Unnecessary parentheses in linear extrapolation. Line 106: `slope = (pb - 14.696)/(1 - 0)` is equivalent to `slope = pb - 14.696`. The `(1-0)` denominator is always 1. Line 107: `intercept = pb - 1*slope` could be `pb - slope`. This is mathematically correct but stylistically unclear; consider simplifying to `slope = pb - 14.696; intercept = 0`.

- **_tables.py:182** — Explicit division by 1000 when Rs column was just created. Line 181-182: Creates `df["Rs (mscf/stb)"]` from rss (scf/stb), then immediately divides by 1000 to convert to mscf/stb. Could combine: `rss_mscf = [r / 1000 for r in rss]` before DataFrame. Minor performance/clarity issue.

- **_density.py:125** — np.log10 without guard on api. Line 125: `np.log10(api)` requires api > 0. While validate_pe_inputs does not explicitly check api, and oil_deno signature requires api OR sg_o, there is a logical path where api=0 is not caught. If both api and sg_o are 0 initially, line 179-180 will raise, but the log10 on line 125 is reached only if sg_sp is undefined. Minor risk; validate_pe_inputs should be enhanced to check api range if supplied.

