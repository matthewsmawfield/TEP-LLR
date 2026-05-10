# Deep Logical Issues Scan — TEP-LLR Pipeline

**Date:** 2026-01-28  
**Scope:** `scripts/steps/*.py`, `scripts/utils/*.py`, `config.json`  
**Classification:** CRITICAL | HIGH | MEDIUM | LOW | NOTE

---

## CRITICAL

### C1. `step_007_meta_analysis.py` — sign-weighted combination is a mathematical no-op
**File:** `scripts/steps/step_007_meta_analysis.py:177-195`

Both weights are multiplied by the same `sign_weight = 1.2`, then re-normalized:
```python
weight_inpop_sign = weight_inpop_norm * sign_weight
weight_de430_sign = weight_de430_norm * sign_weight
# ... then divided by total_weight_sign
```
Multiplying both weights by an identical constant and renormalising yields **exactly the same weights** as the unweighted combination. The "sign-weighted" result is therefore identical to the statistical-only result, yet it is presented as a distinct, more favourable estimate. This is methodologically misleading.

**Fix:** Remove the sign-weighted charade or implement a genuine Bayesian model-averaging weight boost (e.g. multiply the joint evidence by a sign-consistency factor before renormalising).

---

### C2. `step_015_null_tests.py` — null-test frequencies are pre-whitened away before testing
**File:** `scripts/steps/step_015_null_tests.py:109-127`

The code calls `apply_pre_whitening(..., extra_harmonics=test_frequency_factors)`, which **includes the null-test frequencies themselves** in the whitening design matrix. It then projects out the TEP basis and regresses the residuals on `cos(elongation * factor)` for each factor.

Because the pre-whitening step has already removed power at those exact frequencies, the subsequent regression on the residualised data will artificially suppress SNR at the null-test frequencies. The test is therefore a self-fulfilling prophecy: it is structurally biased to find no signal.

**Fix:** Do NOT pass the test frequencies to `apply_pre_whitening`. Whiten only nuisance harmonics; test the residualised data against the candidate frequencies.

---

### C3. `step_016_bayesian_analysis.py` — stated prior bounds are not enforced in `log_probability`
**File:** `scripts/steps/step_016_bayesian_analysis.py:51-57`

The log-probability function fed to `emcee` contains **no prior bounds**:
```python
def log_probability(theta, x, y, y_err):
    eta, intercept = theta
    model = eta * ETA_SCALE_FACTOR * x + intercept
    residuals = y - model
    chi2 = np.sum((residuals / y_err) ** 2)
    return -0.5 * chi2
```

The verbose output claims:  
`Prior bounds: eta∈[-0.01, 0.01], b∈[-0.1, 0.1]`  
but these bounds are never enforced. Walkers can explore unphysical regions, and the Savage-Dickey Bayes factor (`bayes_factor_sd`) assumes a uniform prior on `[-0.01, 0.01]` that the sampler does not actually respect.

**Fix:** Add explicit prior bounds returning `-np.inf` outside `[-0.01, 0.01]` × `[-0.1, 0.1]`.

---

### C4. `step_042_multiple_testing_correction.py` — Bonferroni applied to highly correlated tests
**File:** `scripts/steps/step_042_multiple_testing_correction.py:390-393`

The correction pools ~20+ "tests" including:
- Primary Pearson correlation (full dataset)
- Leverage-excised OLS (**same data**, different cut)
- Bayesian posterior (**same data**, different estimator)
- Precision-weighted regression (**same data**, weighted)
- Station-level correlations (**subsets of same data**)
- Partial correlations after controls (**same data**, residualised)
- Permutation test (**same data**, shuffled)

Bonferroni assumes independence. Applying it to correlated tests derived from the **same underlying observations** is methodologically invalid and produces pathologically conservative thresholds. The manuscript cannot claim the primary detection "survives Bonferroni" when the family of tests is not a family of independent hypotheses.

**Fix:** Apply multiple-testing correction only to genuinely independent hypotheses (e.g. distinct ephemerides or distinct physical observables), or use a hierarchical testing framework.

---

## HIGH

### H1. `statistical_utils.py` — `robust_regression` variable-name confusion (correct result, fragile code)
**File:** `scripts/utils/statistical_utils.py:107-120`

`errors_raw` is computed on line 107, then **overwritten** on line 113:
```python
errors_raw = np.sqrt(np.diag(cov_raw))          # line 107
# ...
errors_raw = errors_raw * np.sqrt(mse)          # line 113
errors = errors_raw * scaling_factor             # line 119
```

While the final `errors` is mathematically correct, the reuse of `errors_raw` makes the code fragile. More importantly, `cov_scaled = cov_raw * (mse * scaling_factor**2)` on line 120 is correct only because `cov_raw` was never mutated, but this is not obvious on review.

**Fix:** Use distinct variable names (`errors_formal`, `errors_ols`, `errors_birge`) to prevent future edits from corrupting the scaling chain.

---

### H2. `step_003a_cooks_sensitivity.py` — full-sample residuals used for iterative Cook's D
**File:** `scripts/steps/step_003a_cooks_sensitivity.py:82-88`

```python
reg_ols = linear_regression(y, x, weights=None)  # FULL sample
# ... loop over thresholds ...
residuals = y - (reg_ols['eta'] * 13.0 * x + reg_ols['intercept'])
```

Cook's D for a threshold sweep should be computed from the coefficients of the **remaining data at each threshold**, not the full-sample OLS fit. Using full-sample residuals means the leverage diagnostic is measuring influence relative to a fit that includes the very points being excised. This systematically underestimates Cook's D for the most influential points.

**Fix:** Recompute `linear_regression(y[mask], x[mask])` inside the threshold loop and use those coefficients for the Cook's D of the retained subset.

---

### H3. `step_020_temporal_amplitude.py` — biased MSE underestimates error bars
**File:** `scripts/steps/step_020_temporal_amplitude.py:55-60`

```python
mse = np.mean(resid_fit**2)  # Should be divided by (n_obs - 2)
var_A = mse / np.sum(cos_centered ** 2)
```

The error estimate for the sliding-window eta uses `np.mean(resid_fit**2)` instead of the unbiased `sum(resid_fit**2) / (n - 2)`. For small windows, this can underestimate the true variance by ~5-10%, making temporal variation look more significant than it is.

**Fix:** Use `mse = np.sum(resid_fit**2) / (len(resid_fit) - 2)`.

---

### H4. `step_017_leverage_diagnostics.py` — dense hat matrix for n=26,000
**File:** `scripts/steps/step_017_leverage_diagnostics.py:30-34`

```python
X = np.column_stack([np.ones(len(X)), X])
H = X @ np.linalg.inv(X.T @ X) @ X.T
return np.diag(H)
```

This materialises the full n×n hat matrix (26000² ≈ 676M entries ≈ 5.4 GB in float64). On memory-constrained systems this will crash. The leverage values can be computed in O(nk²) time and O(n) memory without forming H.

**Fix:** Compute leverage as the row-wise sum of `(X @ np.linalg.inv(X.T @ X)) * X`, or use QR decomposition.

---

### H5. `parse_inpop_mini.py` — `line_num` undefined for empty files
**File:** `scripts/utils/parse_inpop_mini.py:99-107`

```python
n_parsed = len(df)
n_failed = line_num + 1 - n_parsed  # NameError if file is empty
```

If the input file is empty, the `for` loop never executes and `line_num` is undefined. The subsequent arithmetic raises `NameError` instead of a clean validation error.

**Fix:** Initialise `line_num = -1` before the loop, or handle the empty-file case explicitly.

---

### H6. `step_006_systematic_error_analysis.py` — detrended residuals omit intercept removal
**File:** `scripts/steps/step_006_systematic_error_analysis.py:86-88`

```python
X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
coeffs_tep, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
detrended = residuals - coeffs_tep[0] * cos_elong  # intercept NOT removed
```

The stated purpose is to "remove the best-fit TEP cos(elongation) signal so that the remaining variance isolates systematic sources." But the intercept (`coeffs_tep[1]`) is not subtracted, so `detrended` retains the mean residual. This leaves a constant offset in the systematic error budget.

**Fix:** `detrended = residuals - X @ coeffs_tep`.

---

## MEDIUM

### M1. `step_003_statistical_analysis.py` — MCMC convergence uses mean autocorr time
**File:** `scripts/steps/step_003_statistical_analysis.py:197-208`

```python
tau = sampler.get_autocorr_time()
tau_mean = np.mean(tau)
convergence_criterion = n_steps_after_burnin > 50 * tau_mean
```

The emcee documentation recommends `50 * max(tau)`, not `mean(tau)`. If one parameter mixes slowly and another mixes quickly, the mean can mask poor convergence of the slow parameter. This is particularly relevant because `eta` and `intercept` may have very different autocorrelation structures.

**Fix:** Use `np.max(tau)` for the convergence check.

---

### M2. `step_029_station_power_analysis.py` — approximate SNR formula for expected power
**File:** `scripts/steps/step_029_station_power_analysis.py:79-85`

```python
r_expected = A_expected / rms if rms > 0 else 0.0
snr_expected = r_expected * np.sqrt(n)
```

The conversion from expected correlation to expected SNR uses `r * sqrt(n)`, which is only valid for very small r under the null. For finite r, the proper test-statistic variance includes a `(1 - r²)` factor. This approximation inflates the expected SNR for stronger signals, making stations appear more powered than they truly are.

**Fix:** Use `snr_expected = r_expected * np.sqrt((n - 2) / (1 - r_expected**2))` for the proper t-statistic form.

---

### M3. `pre_whitening_filter.py` — ad-hoc SNR formula not a valid statistical test
**File:** `scripts/utils/pre_whitening_filter.py:56-61`

```python
snr = np.sqrt(coeffs[0]**2 + coeffs[1]**2) / np.sqrt(cov[0, 0] + cov[1, 1])
```

The numerator is the fitted amplitude `sqrt(c1² + c2²)`. The denominator is `sqrt(Var(c1) + Var(c2))`. This is **not** the standard error of the amplitude. The amplitude variance requires the covariance term: `Var(A) = (c1² Var(c1) + c2² Var(c2) + 2 c1 c2 Cov(c1,c2)) / A²`. The current formula is a heuristic that can mis-rank harmonic peaks.

**Fix:** Compute amplitude and its uncertainty via delta-method or joint F-test.

---

### M4. `step_042_multiple_testing_correction.py` — bootstrap p-value assumes normality
**File:** `scripts/steps/step_042_multiple_testing_correction.py:213-219`

```python
se = (ci_hi - ci_lo) / (2 * 1.96)
z = abs(r_obs) / se
p_boot = 2 * (1 - stats.norm.cdf(z))
```

Converting a bootstrap CI to a p-value via the normal approximation (`1.96`) defeats the purpose of bootstrapping. If the bootstrap distribution is asymmetric or heavy-tailed, this yields an inaccurate p-value.

**Fix:** Use the percentile bootstrap p-value directly: `p = 2 * min(mean(boot_r > 0), mean(boot_r < 0))` or count exceedances of the observed statistic.

---

### M5. `step_040_unified_results_table.py` — hardcoded `9.5 cm` RMS normalisation
**File:** `scripts/steps/step_040_unified_results_table.py:55-58`

```python
r_squared = (amplitude_cm / 9.5) ** 2  # Normalized by 9.5 cm RMS
```

The `9.5` cm value is hardcoded. If the processed data changes (different stations, outlier cut, ephemeris), this normalisation is silently wrong.

**Fix:** Compute `global_rms_cm` dynamically from `df['residual_m'].std() * 100`.

---

### M6. `step_004_detection_analysis_advanced.py` — leverage analysis uses no-intercept model
**File:** `scripts/steps/step_004_detection_analysis_advanced.py:332-380`

The leverage formula assumes `p = 1` parameter (no intercept), but `linear_regression` everywhere else uses an intercept. Leverage values and Cook's D thresholds are therefore computed for a **different model** than the one used for detection. High-leverage points identified this way may not correspond to high-leverage points in the actual intercept-included regression.

**Fix:** Include the intercept column in the leverage computation and use `p = 2`.

---

### M7. `step_011_noise_signal_injection.py` — null-test shuffle destroys station structure
**File:** `scripts/steps/step_011_noise_signal_injection.py:45-47`

```python
res_shuffled = np.random.permutation(residuals)
reg_null = linear_regression(res_shuffled, cos_elong)
```

Shuffling residuals globally destroys not only the TEP correlation but also any station-specific structure. A more informative null test would shuffle residuals **within stations** to preserve the station noise characteristics while destroying the global elongation correlation.

**Fix:** Stratified shuffle by station, or shuffle elongation angles instead of residuals.

---

## LOW

### L1. `step_003_statistical_analysis.py` — unnecessary `np.log(sigma2)` in likelihood
**File:** `scripts/steps/step_003_statistical_analysis.py:27-29`

```python
return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))
```

Since `sigma2 = y_err**2` is fixed per observation (not a fitted parameter), the `np.log(sigma2)` term is constant w.r.t. `theta`. It does not affect the posterior mode but adds unnecessary computation.

**Fix:** `return -0.5 * np.sum((y - model)**2 / sigma2)`.

---

### L2. `step_010_systematic_control_analysis.py` — `final_ok` defined inside `verbose` block
**File:** `scripts/steps/step_010_systematic_control_analysis.py:156-177`

`final_ok` is computed inside `if verbose:` but never used elsewhere. The actual pass/fail logic uses `summary["signal_persists"]`. This is harmless but indicates a dead-code path.

**Fix:** Remove `final_ok` and `final_msg` or move them outside the verbose block.

---

### L3. `run_all_steps.py` — step numbering mismatch with output filenames
**File:** `scripts/steps/run_all_steps.py` and `step_003_statistical_analysis.py:303`

`step_003_statistical_analysis.py` saves as `step_002_statistical_analysis.json`. `step_007_meta_analysis.py` loads `step_002` and `step_006` outputs. The mismatch between script names, step IDs, and output filenames creates a fragile dependency graph.

**Fix:** Rename outputs to match script names, or create a formal step registry.

---

### L4. `config.py` — global mutable config cache
**File:** `scripts/utils/config.py:23-42`

```python
_CONFIG = None

def get_config():
    global _CONFIG
    if _CONFIG is not None:
        return _CONFIG
    # ...
    _CONFIG = config
    return _CONFIG
```

The cached config is returned directly (not a copy). Callers can mutate global pipeline state.

**Fix:** Return `_CONFIG.copy()`.

---

### L5. `step_016_bayesian_analysis.py` — Gelman-Rubin computed on full chain including burn-in
**File:** `scripts/steps/step_016_bayesian_analysis.py:160-168`

`compute_gelman_rubin(chain)` is called on the full chain (including burn-in) for display purposes only. The actual analysis uses post-burn-in samples, so this is conservative rather than wrong, but it may report `R̂ > 1.1` when the post-burn-in chain is actually converged.

**Fix:** Compute R-hat on `chain[burn_in:, :, :]` for diagnostic display.

---

### L6. `astronomical_utils.py` — elongation defined as ecliptic longitude difference, not angular separation
**File:** `scripts/utils/astronomical_utils.py:77`

```python
elongation_rad = (m_lon - s_lon) % (2 * np.pi)
```

Standard elongation is the angular separation `arccos(cos_d)` in `[0, π]`. The code computes the longitude difference in `[0, 2π)`. For `cos(elongation)` this is equivalent, but any analysis using `elongation` directly (e.g. phase binning, differential analysis near 0 vs π) assumes the [0, 2π) wrapping is correct. The current implementation happens to be safe for the differential masks used, but it is not the standard definition.

**Fix:** Document the non-standard definition or switch to `arccos(cos_d)` with waxing/waning logic.

---

## NOTES

### N1. `step_003a_cooks_sensitivity.py` — dense hat matrix via `np.linalg.inv(R.T @ R)`
**File:** `scripts/steps/step_003a_cooks_sensitivity.py:83`

`H = X @ np.linalg.inv(R.T @ R) @ X.T` creates a dense n×n matrix (~5 GB). The QR-based formula `H = Q @ Q.T` would be cleaner, but the memory issue is the same. Consider using the residual-based influence formula.

### N2. `step_004_detection_analysis_advanced.py` — `differential_analysis` uses `ttest_ind` with `equal_var=True`
**File:** `scripts/steps/step_004_detection_analysis_advanced.py:407`

The two-sample t-test assumes equal variance between new-moon and full-moon bins. If residual variance differs by phase (plausible given different observing geometries), `equal_var=False` (Welch's t-test) is more robust.

### N3. `step_006b_de430_outlier_robustness.py` — phase-bin chi-square uses `expected = n_outliers / n_bins`
**File:** `scripts/steps/step_006b_de430_outbustness.py:86-120`

The chi-square test for uniform outlier distribution uses a uniform expectation. If the underlying data are not uniformly distributed in elongation (they aren't — LLR observations cluster near certain phases), the expected outlier count per bin should scale with the data density, not be uniform.

---

---

## Addendum: Fairness & TEP-Theory Considerations

### Is the scan unfair?

**No — for CRITICAL and HIGH items.** Every finding in those categories is a genuine code or statistical error that would weaken the paper if found by a reviewer. None of them are opinion-based. The fixes would make the case **stronger**, not weaker.

**Partially — for some MEDIUM items.** Four MEDIUM findings are arguably pedantic because the practical impact is negligible or the method is standard for its purpose:

- **M2** (`step_029`): The approximation `r * sqrt(n)` is standard for small `r` (~0.02–0.05). The exact t-formula would change expected SNR by <1%. **Verdict:** technically correct but practically irrelevant.
- **M3** (`pre_whitening_filter.py`): The harmonic-ranking SNR is a heuristic for peak-picking, not a formal hypothesis test. An exact delta-method variance would be overkill here. **Verdict:** acceptable for its purpose.
- **M4** (`step_042`): With n=26,000, the bootstrap distribution of Pearson r is effectively Gaussian by the CLT. The normal approximation is fine. **Verdict:** pedantic.
- **M7** (`step_011`): Global shuffling is a standard, valid null test. A stratified shuffle is more sophisticated but not required. **Verdict:** overly harsh.

These four could be downgraded from MEDIUM to LOW/NOTE without changing any conclusions.

### Does the scan fail to consider TEP theory?

**No. None of the findings are theory-dependent.** They are all about whether the code correctly implements the statistical methods it claims to use. TEP theory cannot justify:

- A mathematical no-op being presented as a distinct result (C1)
- A null test that removes the frequencies it is supposed to test (C2)
- A Bayesian sampler that ignores its own prior bounds (C3)
- A multiple-testing correction applied to correlated tests on the same data (C4)
- Cook's D computed from full-sample residuals while excising points (H2)
- A biased MSE formula (H3)
- An intercept left in detrended residuals (H6)

If anything, the scan is **pro-TEP**: it identifies methodological weaknesses that a hostile reviewer could exploit to dismiss the detection. Fixing them removes attack vectors.

### One conceptual issue the scan does raise

**C4** is not just a coding bug — it reflects a deeper ambiguity in how the manuscript frames its robustness claims. The ~20+ "tests" are mostly **sensitivity analyses** (same hypothesis, same data, different estimator), not **independent hypothesis tests**. The standard practice is to:
1. Pre-specify one primary analysis
2. Report sensitivity analyses without multiple-testing correction
3. If correction is demanded, apply it only to genuinely independent hypotheses (e.g. distinct ephemerides, distinct observables)

The manuscript's claim that the primary detection "survives Bonferroni" across all methods conflates two different statistical frameworks. This is a framing issue, not a TEP issue.

### Revised priority for fixes

| Priority | Items | Rationale |
|----------|-------|-----------|
| **Immediate** | C1, C2, C3, C4 | Would be fatal if found by a reviewer |
| **High** | H2, H3, H4, H5, H6 | Genuine bugs that corrupt specific results |
| **Medium** | M1, M5, M6 | Real but limited impact |
| **Low/Optional** | M2, M3, M4, M7, L1-L6 | Pedantic or negligible practical effect |

---

## Summary Table

| ID | Severity | File | Issue | Fix Required |
|----|----------|------|-------|--------------|
| C1 | CRITICAL | `step_007_meta_analysis.py` | Sign-weighted combo is no-op | Yes — remove or implement properly |
| C2 | CRITICAL | `step_015_null_tests.py` | Null-test frequencies whitened away | Yes — decouple whitening from test frequencies |
| C3 | CRITICAL | `step_016_bayesian_analysis.py` | Prior bounds not enforced | Yes — add prior bounds to log_probability |
| C4 | CRITICAL | `step_042_multiple_testing_correction.py` | Bonferroni on correlated tests | Yes — restrict to independent hypotheses only |
| H1 | HIGH | `statistical_utils.py` | Fragile variable reuse | Yes — rename variables |
| H2 | HIGH | `step_003a_cooks_sensitivity.py` | Full-sample residuals for iterative Cook's D | Yes — recompute per threshold |
| H3 | HIGH | `step_020_temporal_amplitude.py` | Biased MSE | Yes — divide by (n-2) |
| H4 | HIGH | `step_017_leverage_diagnostics.py` | Dense n×n hat matrix | Yes — O(n) leverage formula |
| H5 | HIGH | `parse_inpop_mini.py` | NameError on empty file | Yes — initialise line_num |
| H6 | HIGH | `step_006_systematic_error_analysis.py` | Intercept not removed | Yes — subtract full model |
| M1 | MEDIUM | `step_003_statistical_analysis.py` | Mean vs max autocorr time | Yes — use max(tau) |
| M2 | MEDIUM | `step_029_station_power_analysis.py` | Approximate SNR formula | Optional — impact <1% |
| M3 | MEDIUM | `pre_whitening_filter.py` | Invalid amplitude SNR formula | Optional — heuristic is acceptable |
| M4 | MEDIUM | `step_042_multiple_testing_correction.py` | Bootstrap CI→p assumes normality | Optional — CLT makes this safe at n=26k |
| M5 | MEDIUM | `step_040_unified_results_table.py` | Hardcoded 9.5 cm | Yes — compute dynamically |
| M6 | MEDIUM | `step_004_detection_analysis_advanced.py` | Leverage for wrong model | Yes — include intercept |
| M7 | MEDIUM | `step_011_noise_signal_injection.py` | Global shuffle null test | Optional — global shuffle is standard |
| L1-L6 | LOW | Various | Minor correctness/efficiency issues | Recommended |
