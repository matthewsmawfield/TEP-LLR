## Unified Results Table

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | N | Method | Status |
|-----------|-----------|----------------|-----|---|--------|--------|
| full_systematic_ols | -4.06 | 6.58 | 6.17σ | 25,445 | OLS with full systematic model (cosD + cos2D + sin_m + cos_m + sin_y + cos_y + const) | PRIMARY RESULT |
| ar1_gls_full_model | -4.47 | 9.68 | 4.62σ | 25,445 | Full-model AR(1) GLS with Cochrane-Orcutt on full design matrix + cluster-robust SE (by station) | ROBUSTNESS CHECK - accounts for temporal autocorrelation with systematic controls |
| ar1_gls_cosd_only | -3.29 | 9.37 | 3.51σ | 25,445 | cosD-only AR(1) GLS with Cochrane-Orcutt + cluster-robust SE (by station) | COMPARISON - cosD-only; superseded by full-systematic model |
| full_sample_ols | -3.18 | 6.05 | 5.25σ | 25,445 | cosD-only OLS with 6σ MAD outlier cleaning (step_003) | SECONDARY - cosD-only baseline |
| bayesian_mcmc | -2.86 | 6.50 | 4.40σ | 25,445 | Ensemble MCMC (32 walkers, 3000 steps) | SECONDARY - consistent with primary |
| leverage_excised_ols | -3.31 | 5.85 | 5.65σ | 25,176 | OLS with Cook's Distance excision (threshold: 4/n) | DIAGNOSTIC - confirms leverage inflation |

### Robust Estimands

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | Method | Status |
|-----------|-----------|----------------|-----|--------|--------|
| theil_sen | -2.04 | N/A | N/A | Median of pairwise slopes | ROBUST LOWER BOUND |
| precision_weighted | -3.50 | 6.589113340030604 | 5.3188644467079325 | WLS with 1/σ² station weights | CROSS-STATION VALIDATION |

### Table A: Station-level regression estimates

| Station | N | η (×10⁻⁴) | ση (×10⁻⁴) | η/ση |
|---------|---|-----------|------------|------|
| APO | 2,595 | -2.40 | 0.87 | 2.77 |
| Grasse | 19,390 | -5.40 | 1.09 | 4.97 |
| Haleakala | 737 | 35.54 | 14.49 | 2.45 |
| Matera | 346 | -0.18 | 8.71 | 0.02 |
| McDonald2 | 3,139 | -5.01 | 3.61 | 1.39 |

### Table B: Station-level phase-correlation diagnostics

| Station | r(R, cos D) | p-value | Phase coverage | RMS (cm) |
|---------|-------------|---------|----------------|----------|
| APO | -0.0543 | 5.68e-03 | Good | 3.16 |
| Grasse | -0.0356 | 6.86e-07 | Good | 9.87 |
| Haleakala | 0.0901 | 1.44e-02 | Biased: mean cos(D)=-0.334 | 13.83 |
| Matera | -0.0011 | 9.83e-01 | Biased: mean cos(D)=-0.026 | 6.19 |
| McDonald2 | -0.0248 | 1.65e-01 | Biased: mean cos(D)=-0.316 | 9.55 |