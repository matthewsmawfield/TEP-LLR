## Unified Results Table

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | N | Method | Status |
|-----------|-----------|----------------|-----|---|--------|--------|
| leverage_excised_ols | -3.31 | 5.84 | 5.67σ | 25,177 | OLS with Cook's Distance excision (threshold: 4/n) | PRIMARY RESULT |
| full_sample_ols | -3.17 | 6.04 | 5.25σ | 25,445 | OLS with 6σ MAD outlier cleaning | SECONDARY - inflated by leverage |
| bayesian_mcmc | -3.18 | 6.06 | 5.26σ | 25,445 | Ensemble MCMC (32 walkers, 3000 steps) | SECONDARY - consistent with primary |

### Robust Estimands

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | Method | Status |
|-----------|-----------|----------------|-----|--------|--------|
| theil_sen | -2.04 | N/A | N/A | Median of pairwise slopes | ROBUST LOWER BOUND |
| precision_weighted | -3.50 | 11.299999999999999 | 3.11 | WLS with 1/σ² station weights | CROSS-STATION VALIDATION |

### Station-Level Results

| Station | η (×10⁻⁴) | Error (×10⁻³) | SNR | N | r | p | Powered? |
|---------|-----------|----------------|-----|---|---|---|---------|
| APO | -2.39 | 2.74 | 0.09 | 2,595 | -0.0543 | 5.69e-03 | False |
| Grasse | -5.39 | 1.10 | 0.49 | 19,390 | -0.0357 | 6.82e-07 | False |
| Matera | -0.13 | 14.00 | 0.00 | 346 | -0.0008 | 9.88e-01 | False |
| McDonald2 | -5.00 | 3.77 | 0.13 | 3,139 | -0.0248 | 1.65e-01 | False |
| Haleakala | 35.50 | 10.50 | 0.34 | 737 | 0.0902 | 1.40e-02 | False |