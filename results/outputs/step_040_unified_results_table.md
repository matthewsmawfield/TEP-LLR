## Unified Results Table

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | N | Method | Status |
|-----------|-----------|----------------|-----|---|--------|--------|
| full_systematic_ols | -4.05 | 6.57 | 6.17σ | 25,445 | OLS with full systematic model (cosD + cos2D + sin_m + cos_m + sin_y + cos_y + const) | PRIMARY RESULT |
| ar1_gls_full_model | -4.46 | 9.67 | 4.62σ | 25,445 | Full-model AR(1) GLS with Cochrane-Orcutt on full design matrix + cluster-robust SE (by station) | ROBUSTNESS CHECK - accounts for temporal autocorrelation with systematic controls |
| ar1_gls_cosd_only | -3.28 | 9.36 | 3.51σ | 25,445 | cosD-only AR(1) GLS with Cochrane-Orcutt + cluster-robust SE (by station) | COMPARISON - cosD-only; superseded by full-systematic model |
| full_sample_ols | -3.17 | 6.04 | 5.25σ | 25,445 | cosD-only OLS with 6σ MAD outlier cleaning (step_003) | SECONDARY - cosD-only baseline |
| bayesian_mcmc | -2.87 | 6.61 | 4.35σ | 25,445 | Ensemble MCMC (32 walkers, 3000 steps) | SECONDARY - consistent with primary |
| leverage_excised_ols | -3.31 | 5.84 | 5.67σ | 25,177 | OLS with Cook's Distance excision (threshold: 4/n) | DIAGNOSTIC - confirms leverage inflation |

### Robust Estimands

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | Method | Status |
|-----------|-----------|----------------|-----|--------|--------|
| theil_sen | -2.04 | N/A | N/A | Median of pairwise slopes | ROBUST LOWER BOUND |
| precision_weighted | -3.50 | 11.266104266393777 | 3.1050542288630756 | WLS with 1/σ² station weights | CROSS-STATION VALIDATION |

### Table A: Station-level regression estimates

| Station | N | η (×10⁻⁴) | ση (×10⁻⁴) | η/ση |
|---------|---|-----------|------------|------|
| APO | 2,595 | -2.39 | 0.86 | 2.77 |
| Grasse | 19,390 | -5.39 | 1.09 | 4.97 |
| Haleakala | 737 | 35.48 | 14.46 | 2.45 |
| Matera | 346 | -0.13 | 8.68 | 0.02 |
| McDonald2 | 3,139 | -5.00 | 3.60 | 1.39 |

### Table B: Station-level phase-correlation diagnostics

| Station | r(R, cos D) | p-value | Phase coverage | RMS (cm) |
|---------|-------------|---------|----------------|----------|
| APO | -0.0543 | 5.69e-03 | Biased: mean cos(D)=-0.231 | 3.16 |
| Grasse | -0.0357 | 6.82e-07 | Good | 9.87 |
| Haleakala | 0.0902 | 1.44e-02 | Biased: mean cos(D)=-0.335 | 13.83 |
| Matera | -0.0008 | 9.88e-01 | Biased: mean cos(D)=-0.026 | 6.19 |
| McDonald2 | -0.0248 | 1.65e-01 | Biased: mean cos(D)=-0.316 | 9.55 |