# Correctness — dca

## Critical

### 1. duong_cum: Negative cumulative production for t < 0.001
**Location**: Line 293-295
**Issue**: When `t < 0.001`, `np.linspace(0.001, ti, ...)` generates an array in **reverse order** (descending) because `start > stop`. This causes `np.trapezoid()` to compute a negative area.

**Impact**: High — cumulative production must always be ≥ 0. Function violates physical constraint.

**Example**:
```python
dca.duong_cum(1000, 1.5, 1.2, 0.0001)  # Returns -1.285e-07 (NEGATIVE!)
```

**Root cause**: Hardcoded `0.001` lower bound without handling `ti < 0.001`.

**Fix**: Use `np.linspace(min(0.001, ti), max(0.001, ti), ...)` or handle the `ti < 0.001` case explicitly.

---

### 2. forecast: dt=0 causes ZeroDivisionError
**Location**: Line 925
**Issue**: `t = np.arange(dt, t_end + dt/2, dt)`. When `dt=0`, this raises `ZeroDivisionError` instead of validating input.

**Impact**: Medium — predictable crash without useful error message.

**Fix**: Add `if dt <= 0: raise ValueError("dt must be positive")`

---

### 3. forecast: Accepts uptime outside [0, 1] without warning
**Location**: Line 933
**Issue**: `q = q_capacity * uptime` applies scaling without bounds check. Negative uptime gives negative production; uptime > 1 amplifies beyond capacity.

**Impact**: Medium — produces non-physical results (negative rates) without validation.

**Example**:
```python
fc = dca.forecast(result, t_end=10, uptime=-0.1)  # EUR = -767.4 (nonsensical)
```

**Fix**: Add `if uptime < 0 or uptime > 1: raise ValueError("uptime must be in [0, 1]")`

---

## High

### 4. Hyperbolic fit: Grid bounds miss edge cases (b < 0.05 and b > 0.95)
**Location**: Lines 398, 527
**Issue**: Grid search over `b ∈ [0.05, 0.95]` step 0.01 misses:
  - Very steep hyperbolic decline (b < 0.05, close to exponential)
  - Near-harmonic cases (b > 0.95, close to b=1)

**Impact**: Medium — data with b=0.99 will be fit to b=0.95 instead, with ~0.3% error (acceptable but not optimal).

**Example**:
```python
# True data with b=0.99
result = fit_decline(t, q_true, method='hyperbolic')
# Returns b=0.95 instead of 0.99
```

**Note**: Test `test_fit_hyperbolic_outliers` accepts `0.2 < result.b < 0.9`, which is intentionally loose. The grid is empirically good enough (R² > 0.95) but theoretically incomplete.

---

### 5. EUR: Division by zero when di → 0
**Location**: Line 328
**Issue**: `t_end = -np.log(q_min / qi) / di`. When `di=0`, this raises `ZeroDivisionError`. The code validates `di > 0` at line 321, but the error is not caught gracefully.

**Actually**: Input validation DOES prevent this (di > 0 check), so this is **NOT a bug** (marked as False Alarm).

---

### 6. fit_decline_cum: Non-monotone Np accepted without validation
**Location**: Lines 568-669
**Issue**: Code assumes cumulative production is monotone increasing. If data has shutdowns (Np constant or declining in intervals), the fit may produce invalid results without error.

**Impact**: Low-Medium — edge case (requires bad input data), but silent failure.

**Example**:
```python
Np = np.array([0, 10, 20, 20, 30, 40])  # Flat period
q = np.array([100, 95, 90, 85, 80, 75])
result = fit_decline_cum(Np, q, method='exponential')
# Returns di=0.65, b=0 (incorrect)
```

---

## Medium

### 7. duong_cum: Hardcoded t_min=0.001 ignores time unit
**Location**: Line 293
**Issue**: `t_fine = np.linspace(0.001, ti, ...)` assumes all times > 0.001. This is unit-agnostic (no internal conversions), so 0.001 could mean:
  - Years: 0.36 days (reasonable)
  - Months: 0.03 days (very short)
  - Days: 86 seconds (unrealistic)

**Impact**: Low-Medium — likely not an issue in practice (typical datasets have t >> 0.001 in any unit), but violates unit-agnostic principle.

**Fix**: Remove hardcoded 0.001, or document explicitly.

---

### 8. B-exponent bounds [0, 1] are standard but restrictive
**Location**: Lines 180, 217, 324
**Issue**: Arps equations reject `b > 1`, but literature permits `b ∈ [0, 2]` for unconventionals. Code is intentionally conservative.

**Impact**: Low — affects only niche use cases (ultra-steep decline wells). Conservative choice is defensible.

**Note**: This is a design choice, not a bug. Documented in docstrings.

---

### 9. Duong 'a' parameter: Requires a > 0 but literature permits a < 0
**Location**: Lines 252-253, 283-284, 442-443
**Issue**: Code enforces `a > 0`, but original Duong model can have negative 'a' values. Bounds in `curve_fit` restrict to `[0.01, 10.0]`.

**Impact**: Low — rejecting a < 0 prevents some unconventional fits, but may be intentional simplification.

---

### 10. forecast: Cumulative computed as Riemann sum, not analytical integral
**Location**: Line 948
**Issue**: `Qcum = np.cumsum(q * dt)` uses left-endpoint Riemann approximation. Error is ~2.5% for dt=1 but improves to 0.25% for dt=0.1.

**Impact**: Low-Medium — acceptable for typical use (dt=1 common), but not exact.

**Expected behavior**: Should match `arps_cum()` analytically. Currently requires fine dt for accuracy.

---

## Low

### 11. Constant-rate data: Fit returns di ≈ 0 with R²=0 instead of error
**Location**: Lines 715-737
**Issue**: When data is flat (constant q), exponential fit recovers qi ≈ constant and di ≈ 0, with R²=0. No warning that decline signal is absent.

**Impact**: Very Low — edge case, R²=0 signals poor fit, but silent acceptance could mask input errors.

**Example**:
```python
q = [100, 100, 100, 100]
result = fit_decline(t, q, method='exponential')
# Returns di=0.00001, R²=0.0 (flat, useless)
```

---

### 12. fit_decline window: Shifted time may confuse users
**Location**: Lines 696-706
**Issue**: When `t_start` is set, returned `qi` is the rate at the window start time, not at t=0. This is mathematically correct but unintuitive.

**Impact**: Very Low — documented in docstring ("Shift so window starts at t=0").

---

### 13. ratio_forecast: No validation of domain
**Location**: Lines 868-896
**Issue**: `ratio_forecast()` calls itself recursively via `forecast()` without checking that `x` (time or cumulative) is compatible with `result.domain`. Mismatches silently evaluate at wrong x-axis.

**Actually**: Not a bug—forecast() handles domain routing correctly. *False Alarm.*

---

### 14. Forecast with q_min very close to initial rate
**Location**: Lines 936-946
**Issue**: If `q_min > q[0]` (economic limit above initial rate), the array is truncated to first point only: `t = t[:1]`, giving 1-point forecast with eur ≈ q[0] * dt.

**Impact**: Very Low — edge case, mathematically correct (no production above limit), but may surprise users.

---

## Summary

| Severity | Count | Issue |
|----------|-------|-------|
| Critical | 3 | duong_cum negative, dt=0 crash, uptime unbounded |
| High | 2 | Hyperbolic grid bounds, Duong a > 0 constraint |
| Medium | 4 | duong_cum t_min unit-agnostic, b bounds restrictive, non-monotone Np, Duong a < 0 |
| Low | 4 | Constant q returns flat result, window shifting opacity, q_min edge case, etc. |

**Regression risk**: Very Low. All 63 existing tests pass. Issues are edge cases or design constraints.

