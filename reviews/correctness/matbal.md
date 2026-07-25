# Correctness — matbal

## Critical

### 1. Gas P/Z Metric pvt_table Interpolation Bug
**Location**: Lines 250, 256  
**Severity**: Critical — incorrect results in metric mode with tabulated PVT  
**Issue**: When `metric=True` and `pvt_table` provided, code converts p and pvt_p to psia (lines 235-236) but then interpolates using the original p and pvt_p (in bars):
```python
if metric:
    p_psia = p * BAR_TO_PSI  # p_psia created
    pvt_p_psia = pvt_p * BAR_TO_PSI
else:
    p_psia = p
    pvt_p_psia = pvt_p
# ...later:
z = np.interp(p, pvt_p, pvt_Z)  # BUG: should use p_psia, pvt_p_psia
bg = np.interp(p, pvt_p, pvt_Bg)  # BUG: should use p_psia, pvt_p_psia
```
When metric=True, p and pvt_p are in bars but the table pressures (from line 232 onwards) are in psia for computation. Interpolating bars onto a psia-scaled table gives nonsense.

**Fix**: Replace lines 250 and 256 with:
```python
z = np.interp(p_psia, pvt_p_psia, pvt_Z)
bg = np.interp(p_psia, pvt_p_psia, pvt_Bg)
```

**Test coverage**: `test_gas_matbal_pvt_z` and `test_gas_matbal_pvt_bg` only test in oilfield units (metric=False). No metric pvt_table tests exist.

---

## High

### 2. Oil Havlena-Odeh Efw Multiplier — (1+m) vs Standard Form
**Location**: Line 551, 604, 642  
**Severity**: High — may violate standard material balance formulation  
**Issue**: The code uses:
```python
denom = Eo + m_t * Eg + (1.0 + m_t) * Efw
```
Standard Havlena-Odeh from Dake/Craft-Hawkins is:
```
F = N*Eo + N*m*Eg + N*Efw
```
not `N*(1+m)*Efw`. The (1+m) multiplier suggests either:
- Efw is being double-counted (once for oil zone, once for gas cap), OR
- This is a non-standard formulation variant that should be documented

**Impact**: Changes OOIP estimate when formation/water compressibility significant.  
**Verification needed**: Compare against published Havlena-Odeh examples from SPE papers (e.g., Havlena-Odeh 1963, Craft et al. "Applied Petroleum Reservoir Engineering").

**Workaround**: Document explicitly if intentional; otherwise, remove the (1+m) and use just Efw.

---

### 3. Regression Objective — Silent Fallback Failure Risk
**Location**: Lines 586-595  
**Severity**: High — Rust acceleration failure mode  
**Issue**: When Rust available and called, if it raises ImportError or AttributeError, execution continues to Python fallback. However, the try/except swallows the exception without logging:
```python
if _RUST_AVAILABLE:
    try:
        return _rust.matbal_objective_rust(...)  # Returns early if success
    except (ImportError, AttributeError):
        pass  # Silent swallow — no warning logged
# Falls through to Python version
```
If the Rust library is corrupted or partially loaded, the fallback is silent and users don't know they're running Python instead of accelerated code.

**Fix**: Add explicit fallthrough or warning:
```python
try:
    return _rust.matbal_objective_rust(...)
except (ImportError, AttributeError) as e:
    import warnings
    warnings.warn(f"Rust matbal_objective unavailable, using Python fallback: {e}")
```

---

### 4. No Validation of Optional Rp Parameter
**Location**: Lines 471-476  
**Severity**: High — silent error if Rp wrong length  
**Issue**: When `Rp` provided, code assumes it matches p length but doesn't validate:
```python
if Rp is None:
    Rp_arr = None
else:
    Rp_arr = np.asarray(Rp, dtype=float)
    if metric:
        Rp_arr = Rp_arr * SM3_PER_SM3_TO_SCF_PER_STB
# No length check!
```
If Rp has wrong length, error occurs later during array operations (e.g., line 533: `Rp_arr - Rs_arr`).

**Fix**: Add after line 474:
```python
if len(Rp_arr) != n:
    raise ValueError(f"Rp must have same length as p, got {len(Rp_arr)} and {n}")
```

---

## Medium

### 5. Drive Indices Not Guaranteed to Sum to 1.0
**Location**: Lines 637-642  
**Severity**: Medium — diagnostic accuracy  
**Issue**: DDI, SDI, CDI are computed as:
```python
ddi = N_pts * Eo / F_safe
sdi = N_pts * m * Eg / F_safe
cdi = N_pts * (1.0 + m) * Efw / F_safe
```
Test `test_oil_matbal_drive_indices_sum` checks they sum to ~1.0, but only within tolerance `< 0.05`. The (1+m) multiplier on CDI means:
```
DDI + SDI + CDI = (Eo + m*Eg + (1+m)*Efw) / denom ≠ 1.0 always
```
If the formulation is indeed non-standard (see Issue #2), drive indices won't represent true fractional contribution to drive.

**Impact**: Users may misinterpret drive mechanism fractions.  
**Workaround**: Document that indices may not sum to 1.0 due to formation/water compressibility accounting.

---

### 6. Cole Plot (F/Et) Interpretation Assumes Volumetric
**Location**: Lines 299-302, test at line 108-121  
**Severity**: Medium — diagnostic misuse risk  
**Issue**: Code computes `cole_F_over_Et = F / Et` and test checks it's approximately constant. This is only valid for **volumetric reservoirs** (no aquifer, no compressibility). If aquifer or compaction present, F/Et varies even if reservoir is depletion-drive.

**Test gap**: `test_gas_matbal_cole_volumetric` is named to indicate it's for volumetric case, but code comment says "Cole plot diagnostic" without caveat. Users unfamiliar with Cole plot methodology may misuse.

**Documentation**: Add warning in docstring and test that Cole plot assumes volumetric behavior.

---

## Low

### 7. Regression Bounds Conversion — Clarity Issue
**Location**: Lines 568-572  
**Severity**: Low — confusion risk but working  
**Issue**: When `metric=True` and `regress` provided with cf/cw bounds, code converts bounds from 1/bar to 1/psi for internal optimization, then converts result back. This is correct but non-obvious:
```python
if metric:
    bounds = [
        (lo / BAR_TO_PSI, hi / BAR_TO_PSI) if k in ('cw', 'cf') else (lo, hi)
        for k, (lo, hi) in zip(param_names, bounds)
    ]
```
Later at line 627: `regressed_result[k] = float(v * BAR_TO_PSI)`. Easy to get backwards.

**Suggestion**: Add comment: `# Convert from user's 1/bar to internal 1/psi units`

---

### 8. Bg Calculation in Oil Matbal Not Tagged as rb/scf
**Location**: Line 520  
**Severity**: Low — unit clarity  
**Issue**: When computing Bg for oil_matbal, code does:
```python
bg_i = gas.gas_bg(pi, sg_sp, degf_field) / CUFTperBBL  # comment: converts to rb/scf
```
Bg from gas module is in rcf/scf; dividing by 5.61458 gives rb/scf. This is correct but unusual. Standard material balance uses Bg in rcf/scf; the conversion to rb/scf is to match F units (which are in rb).

**Suggestion**: Comment should clarify: `# Convert rcf/scf → rb/scf for consistency with F (in rb)`

---

### 9. No Upper Bound Check on Gp/Np
**Location**: Input validation (lines 195-203, 405-414)  
**Severity**: Low — edge case  
**Issue**: Code validates pressures are positive and arrays match length, but doesn't check if production is monotonically increasing or bounded. Negative or decreasing Gp/Np could indicate data entry error but is silently accepted.

**Current behavior**: Code will compute; if Gp or Np decreases, denom or numerator may become negative, giving nonsense results.

**Suggestion**: Add optional check (non-blocking, since some workflows may use backward production):
```python
if np.any(np.diff(Np) < -1e-6):  # Allow small numerical noise
    warnings.warn("Np not monotonically increasing — possible data error")
```

---

## Summary

**Counts**: 1 Critical, 4 High, 2 Medium, 3 Low

**Top 3 Blockers**:
1. **Gas pvt_table metric interpolation bug** — Wrong results in metric mode; test coverage gap.
2. **Havlena-Odeh (1+m)*Efw formulation** — Possible deviation from standard; needs verification against reference papers.
3. **Silent Rust fallback failure** — Users unaware of performance degradation; add diagnostic warning.

All other issues are documentation, edge-case handling, or assumptions that should be made explicit.
