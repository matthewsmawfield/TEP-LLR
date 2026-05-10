# TEP-LLR Full Codebase Audit Plan

## Document Purpose

Records findings of the full codebase audit. Each issue: root cause, impact on scientific inference, concrete fix, and priority.

---

## Section 1: Previously Fixed Issues — Verified Correct

| ID | File | Fix | Status |
|---|---|---|---|
| C1 | `step_007_meta_analysis.py` | Sign-weighted no-op removed; sign consistency qualitative only | Verified |
| C2 | `step_015_null_tests.py` | Null-test frequencies decoupled from pre-whitening | Verified |
| C3 | `step_016_bayesian_analysis.py` | Prior bounds `[-0.01, 0.01]` enforced in `log_probability` | Verified |
| C4 | `step_042_multiple_testing_correction.py` | Distinguishes 4 independent hypotheses from 14 sensitivity analyses | Verified |
| H1 | `statistical_utils.py` | Variables renamed (`errors_formal`, `mse_unbiased`, `errors_ols`, `errors_birge`) | Verified |
| H2 | `step_003a_cooks_sensitivity.py` | Regression recomputed per threshold on masked data | Verified |
| H3 | `step_020_temporal_amplitude.py` | Unbiased MSE: `sum(resid^2)/(n-2)` | Verified |
| H4 | `step_017_leverage_diagnostics.py` | O(n) memory leverage computation | Verified |
| H5 | `parse_inpop_mini.py` | `line_num = -1` initialised before loop | Verified |
| H6 | `step_006_systematic_error_analysis.py` | Full model (slope + intercept) subtracted for detrending | Verified |

---

## Section 2: New Critical Issues

### C5: Dense Hat Matrix in `step_025_leverage_temporal_clustering.py`

**Location:** `scripts/steps/step_025_leverage_temporal_clustering.py:34`

**Root Cause:**
```python
H = X @ np.linalg.inv(X.T @ X) @ X.T
leverage = np.diag(H)
```
This materialises the full n x n hat matrix. For n = 26,000, H has ~676M entries (~5.4 GB float64), crashing memory-constrained systems.

**Impact:** Step 025 aborts; temporal clustering diagnostic is missing from the audit trail.

**Fix:** Use O(n) computation (same as H4 fix in step_017):
```python
XtX_inv = np.linalg.inv(X.T @ X)
leverage = np.sum((X @ XtX_inv) * X, axis=1)
```

**Priority: 1 (Immediate)**

---

## Section 3: New High-Severity Issues

### H7: Leverage Computed for Inconsistent Model in `step_004_detection_analysis_advanced.py`

**Location:** `scripts/steps/step_004_detection_analysis_advanced.py:332-338`

**Root Cause:**
```python
# For single predictor WITHOUT intercept
x_sq_sum = np.sum(cos_elong**2)
leverage = cos_elong**2 / x_sq_sum
```
The comment explicitly says "without intercept", but `linear_regression` always adds an intercept column. Leverage and Cook's D are computed for a regression-through-origin model, inconsistent with the actual fitted model. Points near cos_elong = 0 get zero leverage regardless of their influence on the intercept term.

**Impact:** Incorrect influence diagnostics; downstream robustness claims may be based on flawed thresholds.

**Fix:** Use the same design matrix as `linear_regression`:
```python
X = np.column_stack([np.ones(n), cos_elong])
XtX_inv = np.linalg.inv(X.T @ X)
leverage = np.sum((X @ XtX_inv) * X, axis=1)
```

**Priority: 2**

---

### H8: Biased MSE in IPW Error Computation in `step_021_ipw_validation.py`

**Location:** `scripts/steps/step_021_ipw_validation.py:102-103`

**Root Cause:**
```python
resid_w = yw - Xw @ coeffs
mse = np.sum(resid_w**2) / len(resid_w)
```
For a 2-parameter weighted model, MSE should divide by effective degrees of freedom (not raw n). Dividing by `len(resid_w)` underestimates variance, inflating IPW significance.

**Impact:** IPW test may pass spuriously; key manuscript robustness argument is weakened.

**Fix:** Use effective sample size for weighted regression:
```python
n_eff = np.sum(weights)**2 / np.sum(weights**2)
dof = n_eff - 2
mse = np.sum(resid_w**2) / dof if dof > 0 else np.nan
```

**Priority: 3**

---

### H9: Crude Solar Cycle Approximation in `step_023_solar_cycle_correlation.py`

**Location:** `scripts/steps/step_023_solar_cycle_correlation.py:26`

**Root Cause:**
```python
def solar_activity_index(years_array):
    return 0.5 * (1.0 + np.cos(2 * np.pi * (years_array - 1990.0) / 11.0))
```
Real solar activity (sunspot number) has variable period (9.5-13 yr), asymmetric rise/fall, and amplitude modulation. A rigid cosine anchored to 1990 is systematically misaligned with actual solar extrema during 1985-2020.

**Impact:** High/low solar bins contain epochs that are not actually at solar extrema; correlation results are unreliable.

**Fix:** Replace with actual SILSO sunspot number data interpolated to observation epochs, or a spline fit to known minima/maxima dates.

**Priority: 4**

---

### H10: Inconsistent Prior Bound Style in `step_003_statistical_analysis.py`

**Location:** `scripts/steps/step_003_statistical_analysis.py:33-37`

**Root Cause:**
```python
if -0.01 < eta < 0.01 and -0.1 < intercept < 0.1:
```
Uses strict inequalities (`<`). `step_016_bayesian_analysis.py:58` uses inclusive bounds (`<=`). Inconsistency between two Bayesian analyses that should share the same prior.

**Impact:** Low direct inference impact, but damages reproducibility and methodological rigour.

**Fix:** Change to inclusive bounds:
```python
if -0.01 <= eta <= 0.01 and -0.1 <= intercept <= 0.1:
```

**Priority: 5 (one-line fix, high rigour payoff)**

---

### H11: Nonsensical Fallback Correlation and Latent NameError in `step_007_meta_analysis.py`

**Location:** `scripts/steps/step_007_meta_analysis.py:92-96`

**Root Cause:**
```python
de430_r = de430_data['eta'] / de430_data['eta_error'] * \
          (de430_data['eta_error'] / 0.266)
```
Simplifies to `eta / 0.266` (dimensionless / metres = physically meaningless). Also, `de430_residuals` is undefined in the `else` branch; `'de430_residuals' in locals()` is always False, so `de430_p` is always `None`.

**Impact:** Latent crash or incomplete JSON if downstream code does not handle `None`.

**Fix:** Remove fallback. Fail loudly if upstream data is missing required fields:
```python
if 'eta' not in de430_data or 'eta_error' not in de430_data:
    raise ValueError("step_006 DE430 comparison missing required fields.")
```

**Priority: 6**

---

### H12: Bootstrap p-Value Assumes Normality in `step_042_multiple_testing_correction.py`

**Location:** `scripts/steps/step_042_multiple_testing_correction.py:213-216`

**Root Cause:**
```python
se = (ci_hi - ci_lo) / (2 * 1.96)
z = abs(r_obs) / se
p_boot = 2 * (1 - stats.norm.cdf(z))
```
Back-calculates SE from percentile span assuming a normal bootstrap distribution. Bootstrap distributions of correlations are often asymmetric or heavy-tailed. Normal approximation can be anti-conservative.

**Impact:** Bootstrap p-value may be unreliable; could push marginal tests below threshold spuriously.

**Fix:** Use non-parametric p-value from bootstrap sample directly:
```python
n_less = np.sum(boot_r < 0)
n_greater = np.sum(boot_r > 0)
p_boot = 2 * min(n_less, n_greater) / len(boot_r)
```

**Priority: 7**

---

## Section 4: New Medium-Severity Issues

### M1: Misleading "Ephemeris-Independent" Label in `step_009`

**Location:** `scripts/steps/step_009_ephemeris_independent_analysis.py`

**Root Cause:** Title claims "ephemeris-independent" but the step still uses `elongation_rad` computed from Skyfield/DE421 in preprocessing. There is no independent elongation calculation.

**Impact:** Misleading label damages credibility if reviewers inspect the code.

**Fix:** Rename to "Systematic-Corrected Residual Analysis" and update docstring.

**Priority: 8**

---

### M2: Inconsistent Regression Tool in `step_027_day_night_thermal_bias.py`

**Location:** `scripts/steps/step_027_day_night_thermal_bias.py:104,111,117`

**Root Cause:** Uses `statsmodels.api.OLS` directly instead of the project's `linear_regression` utility, which applies Birge scaling, condition-number checks, and unbiased MSE.

**Impact:** Slight numerical inconsistency; harder to maintain.

**Fix:** Refactor to use `linear_regression` from `scripts.utils.statistical_utils`.

**Priority: 9**

---

### M3: AR(1) Rho Error Assumes White Noise in `step_003_statistical_analysis.py`

**Location:** `scripts/steps/step_003_statistical_analysis.py:77-81`

**Root Cause:** Uses Bartlett formula `sqrt((1 - rho^2) / (n - 2))` for AR(1) coefficient SE. This formula is only valid for white-noise residuals. The proper asymptotic SE for AR(1) is `sqrt((1 - rho^2) / n)`.

**Impact:** Underestimates uncertainty in rho; may mislead about autocorrelation significance.

**Fix:** Replace with asymptotic SE: `se_rho = np.sqrt((1 - rho**2) / n)`.

**Priority: 10**

---

### M4: Docstring / Filename Mismatch in `step_006_multi_ephemeris_comparison.py`

**Location:** `scripts/steps/step_006_multi_ephemeris_comparison.py:3`

**Root Cause:** Docstring says "Step 005" but file is `step_006` and pipeline lists it as step 006.

**Fix:** Update docstring to "Step 006".

**Priority: 11**

---

### M5: Non-Standard Effective Sample Size in `step_019_systematic_monte_carlo.py`

**Location:** `scripts/steps/step_019_systematic_monte_carlo.py:163`

**Root Cause:**
```python
n_eff = np.sum(cos_centered**2) / np.max(cos_centered**2)
```
Denominator depends on the single most extreme observation, making `n_eff` unstable.

**Impact:** Unstable `eta_error` weights can distort the systematic error distribution.

**Fix:** Use standard formula or simply use `n` (negligible difference for n=26,000), or delegate to `linear_regression`.

**Priority: 12**

---

## Section 5: Recommended Implementation Order

| Priority | Issue | Rationale |
|---|---|---|
| 1 | C5 | Prevents pipeline crash; memory safety is foundational. |
| 2 | H7 | Incorrect diagnostics undermine robustness narrative. |
| 3 | H8 | Inflates significance of key IPW test. |
| 4 | H9 | Wrong solar proxy makes correlation results unreliable. |
| 5 | H10 | One-line fix; removes prior inconsistency. |
| 6 | H11 | Latent crash bug if upstream data changes. |
| 7 | H12 | Anti-conservative p-values weaken multiple-testing correction. |
| 8 | M1 | Misleading label damages credibility. |
| 9 | M2 | Numerical inconsistency; should use shared utilities. |
| 10 | M3 | Underestimates AR(1) uncertainty. |
| 11 | M4 | Documentation hygiene. |
| 12 | M5 | Unstable error estimate in MC error budget. |

---

## Section 6: Files to Modify

1. `scripts/steps/step_025_leverage_temporal_clustering.py` — C5
2. `scripts/steps/step_004_detection_analysis_advanced.py` — H7
3. `scripts/steps/step_021_ipw_validation.py` — H8
4. `scripts/steps/step_023_solar_cycle_correlation.py` — H9
5. `scripts/steps/step_003_statistical_analysis.py` — H10, M3
6. `scripts/steps/step_007_meta_analysis.py` — H11
7. `scripts/steps/step_042_multiple_testing_correction.py` — H12
8. `scripts/steps/step_009_ephemeris_independent_analysis.py` — M1
9. `scripts/steps/step_027_day_night_thermal_bias.py` — M2
10. `scripts/steps/step_006_multi_ephemeris_comparison.py` — M4
11. `scripts/steps/step_019_systematic_monte_carlo.py` — M5
