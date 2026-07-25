# Correctness — recommend

## Critical
None.

## High
None.

## Medium

### 1. `well_type` Parameter Unused in `recommend_vlp_method()`
**Location**: `recommend.py:198-232` (function signature and body)

**Issue**: The `well_type` parameter is declared in both `recommend_vlp_method()` (line 198) and propagated through `recommend_methods()` (line 273), but is never actually used in the function body. The documentation states it should be 'gas' or 'oil', yet the same VLP recommendations are returned regardless of fluid type.

**Impact**: 
- If oil wells require different VLP correlations than gas wells, this parameter could silently ignore that requirement.
- No immediate correctness failure (works today), but masks incomplete implementation.

**Recommendation**: Either (a) remove the unused parameter to clarify intent, or (b) implement fluid-type-aware logic if oil-specific VLP recommendations are needed.

---

## Low

### 1. No Input Validation on Composition Fractions
**Location**: `recommend_gas_methods()` (line 71-151)

**Issue**: Function accepts `co2`, `h2s`, `n2`, `h2` parameters but performs no validation. Does not check:
- Negative fractions (e.g., `co2=-0.1`)
- Sum exceeding 1.0 (e.g., `co2=0.6, h2s=0.6`)

**Current Behavior**: Invalid inputs are silently accepted and may produce incorrect recommendations based on nonsensical composition.

**Rationale**: The `validate_pe_inputs()` function in `shared_fns.py` (lines 176-200) exists specifically for this validation and IS called by `gas_rate_radial()`, `gas_rate_linear()`, and other gas module functions. It enforces:
```python
if frac_sum > 1.0:
    raise ValueError(f"Sum of non-hydrocarbon mole fractions ({frac_sum}) exceeds 1.0")
```

**Recommendation**: Call `validate_pe_inputs()` at the start of `recommend_gas_methods()` with composition parameters. This would catch data entry errors before generating recommendations.

---

### 2. VELAR Correlation Applicability Range Not Documented or Checked
**Location**: `recommend_oil_methods()` (line 154-195) and oil module usage

**Issue**: Function always recommends VELAR for any API value without checking applicability. In contrast, the oil module itself documents validity ranges:
- Standing: API [16.5, 63.8], T [100–258°F]
- VALMC: API [6, 56.8], T [78–330°F]
- VELAR: **No range documented** (despite being the recommended default)

**Current Behavior**: VELAR is recommended as the primary method for all API values (< 10, 10-50, > 50), with identical recommendation regardless of API range.

**Impact**: No immediate failure (VELAR likely extrapolates outside calibration range rather than erroring). However, produces warnings in the calculation layer (`oil.py`) if invoked at extreme API values, creating asymmetry: recommend module is unaware of validity limits.

**Recommendation**: 
1. Investigate and document VELAR's actual applicability range from the original McCain/Blasingame papers.
2. Add a warning in `recommend_oil_methods()` if API falls outside VELAR's documented range, similar to the gas module's Tr/Ppr range warnings.
3. Consider recommending alternative methods (e.g., VALMC for very heavy oils) if VELAR is confirmed to have narrow applicability.

---

## Code Quality Notes (Out of Scope)

These are mentioned for context but are NOT correctness issues:

- Well-typed dataclass for `MethodRecommendation` is clean and reusable.
- Decision tree logic is sound: H2 > high inerts (55%) > moderate CO2/H2S (10%) > clean gas.
- VLP deviation boundary at 30° is correct per industry guidance (HB/Gray are vertical-development correlations).
- All method strings match valid Enum members in `classes.py`.
- All 13 unit tests pass; all doc example tests pass.

