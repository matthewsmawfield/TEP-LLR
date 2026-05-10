## Unified Results Table

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | N | Method | Status |
|-----------|-----------|----------------|-----|---|--------|--------|
| leverage_excised_ols | -3.31 | 5.84 | 5.67σ | 25,177 | OLS with Cook's Distance excision (threshold: 4/n) | PRIMARY RESULT |
| full_sample_ols | -3.17 | 6.04 | 5.25σ | 25,445 | OLS with 6σ MAD outlier cleaning | SECONDARY - inflated by leverage |
| bayesian_mcmc | -2.87 | 6.61 | 4.35σ | 25,445 | Ensemble MCMC (32 walkers, 3000 steps) | SECONDARY - consistent with primary |

### Robust Estimands

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | Method | Status |
|-----------|-----------|----------------|-----|--------|--------|
| theil_sen | -2.04 | N/A | N/A | Median of pairwise slopes | ROBUST LOWER BOUND |
| precision_weighted | -3.50 | 11.266104266393777 | 3.1050542288630756 | WLS with 1/σ² station weights | CROSS-STATION VALIDATION |

### Station-Level Results

| Station | η (×10⁻⁴) | Error (×10⁻³) | SNR | N | r | p | Powered? |
|---------|-----------|----------------|-----|---|---|---|---------|
| APO | -2.39 | 0.09 | 2.77 | 2,595 | -0.0543 | 5.69e-03 | False |
| Grasse | -5.39 | 0.11 | 4.97 | 19,390 | -0.0357 | 6.82e-07 | True |
| Haleakala | 35.48 | 1.45 | 2.45 | 737 | 0.0902 | 1.44e-02 | False |
| Matera | -0.13 | 0.87 | 0.02 | 346 | -0.0008 | 9.88e-01 | False |
| McDonald2 | -5.00 | 0.36 | 1.39 | 3,139 | -0.0248 | 1.65e-01 | False |