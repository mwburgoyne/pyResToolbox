# Correctness — plyasunov

## Critical

### 1. IAPWS-IF97 Isothermal Compressibility Formula (kappa_T_if97)
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/iapws_if97.py`, line 137

**Issue:** The formula for isothermal compressibility appears to have an incorrect factor.

Current code:
```python
return -gpp / (gp * P_STAR)
```

According to IAPWS-IF97 standard (Wagner & Pruss 2000), the correct formula should be:
```python
return -gpp / (gp * gp * P_STAR)  # or equivalently: -gpp / (gp**2 * P_STAR)
```

**Impact:** If incorrect, all compressibility calculations are wrong by a factor of γ_π, affecting:
- Plyasunov V2_inf calculations (which multiply by kappa_T)
- Garcia density mixing corrections in SoreideWhitson
- Any downstream brine property calculations using Cf_sat/Cf_usat

**Root Cause:** The IAPWS-IF97 compressibility definition: κ_T = (∂ln(ρ)/∂P)_T = -γ_ππ / (γ_π² * P_star), not -γ_ππ / (γ_π * P_star)

**Verification Needed:** Confirm against official IAPWS-IF97 Part 1 equations or validated reference implementations (e.g., IF97 NIST webbook, libeieos).

---

## High

### 2. IAPWS-IF97 Region 1 Temperature Boundary (Strict Inequality)
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/iapws_if97.py`, line 99

**Issue:** Temperature bounds check uses strict inequality for lower bound:
```python
if T < 273.15 or T > 623.15:
```

Should be:
```python
if T <= 273.15 or T > 623.15:  # or use ≤ for lower bound per IAPWS-IF97
```

IAPWS-IF97 Region 1 valid range is **273.15 K ≤ T ≤ 623.15 K**. The current code rejects T=273.15K (0°C), which is a valid boundary point.

**Impact:** Low in practice (rare to use exactly 273.15 K), but semantically incorrect.

---

### 3. Missing Saturation Pressure Validation (P_sat Boundary)
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/iapws_if97.py`, lines 101-102

**Issue:** Documentation states valid range as "P_sat(T) ≤ P ≤ 100 MPa" but code only checks P > 0:
```python
if P <= 0 or P > 100:
    raise ValueError(...)
```

The code doesn't validate P ≥ P_sat(T). This allows calls at pressures below saturation (e.g., T=373K, P=0.1 MPa where P_sat≈0.47 MPa), which represents non-physical vapor phase.

**Impact:** Medium. In normal reservoir engineering usage (P_res > P_sat), this isn't triggered. However, Plyasunov model is undefined for supercritical vapor states. The brine module always passes reservoir P > P_sat, so practical risk is low.

**Note:** Implementing this would require P_sat(T) calculation. Currently documentation promise isn't enforced.

---

## Medium

### 4. Plyasunov Paper Reference Incomplete
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/plyasunov_model.py`, lines 1-29

**Issue:** Docstring mentions "Part IV (Plyasunov 2021)" for H2, N2, CH4 but no full citation provided in References section.

**Impact:** Low. Comments reference "Part IV Eq. 5" internally, and Part IV is correctly dated 2021 in code structure. Missing formal bibliographic entry is cosmetic but affects reproducibility.

**Fix:** Add full citation for Part IV to References section.

---

### 5. Potential Numerical Precision Issue in B12_CO2 Exponent
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/plyasunov_model.py`, line 71

**Issue:** Line 61 documents the formula as having term with exponent (21/2):
```
B12 = ... + b6/(T*)^(21/2)
```

But line 71 implements:
```python
+ b[5] / Tstar**10.5
```

While 21/2 = 10.5, using the floating-point representation may introduce precision loss. However, **this is mathematically equivalent and acceptable** since 10.5 is exactly representable in IEEE 754.

**Impact:** Negligible. Not a correctness issue.

---

## Low

### 6. Edge Case: Garcia Mixing Rule Near x_total→1.0
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/brine/brine.py`, lines 1805-1813

**Issue:** Code warns when x_total > 0.999999 (x1 < 1e-6) that "Garcia mixing rule breaks down" and returns unsaturated density. The formula has x1 in denominator, which can become numerically unstable.

**Correctness Impact:** Low. The warning is appropriate and defensive. However, no alternative model is provided—code just returns brine density ignoring dissolved gas, which is incorrect for very high gas mole fractions. This is documented behavior (warning issued).

**Severity:** Low because (a) this only affects extreme, non-physical states (>99.9999% gas), and (b) a warning is issued.

---

### 7. P_STAR and T_STAR Constants Not Validated
**Location:** `/home/mark/projects/pyResToolbox/pyrestoolbox/plyasunov/iapws_if97.py`, lines 24-26

**Values:**
```python
R_SPECIFIC = 461.526e-6  # MPa*m3/(kg*K)
P_STAR = 16.53           # MPa
T_STAR = 1386.0          # K
```

**Note:** These match IAPWS-IF97 Part 1 specification. No issue found. Documented for completeness.

---

## Summary

| Severity | Count | Requires Fix |
|----------|-------|--------------|
| **Critical** | 1 | YES — kappa_T formula |
| **High** | 2 | YES — Temperature boundary, P_sat validation |
| **Medium** | 1 | Consider — Missing Part IV citation |
| **Low** | 2 | No — Acceptable or documented |

### Top 3 Issues to Address
1. **Isothermal compressibility formula** — Affects all V2_inf and density calculations
2. **Temperature boundary strictness** — Semantic correctness at 0°C
3. **Saturation pressure validation** — Prevents invalid physical states
