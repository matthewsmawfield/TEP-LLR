## Unified Results Table

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | N | Method | Status |
|-----------|-----------|----------------|-----|---|--------|--------|
| precision_weighted_full_systematic | -3.91 | 5.63 | 6.94σ | 25,445 | Full-systematic WLS with 1/σ² per-station weights on cosD + cos2D + sin_m + cos_m + sin_y + cos_y + const | PRIMARY HEADLINE ESTIMAND — uses every observation; no data deleted |
| cooks_excised_full_systematic | -3.87 | 4.95 | 7.82σ | 23,837 | Full-systematic OLS with Cook's Distance excision (D > 4/n) on cosD + cos2D + sin_m + cos_m + sin_y + cos_y + const | SECONDARY ROBUSTNESS CHECK — confirms signal persists after removing 1,608 high-leverage points |
| full_systematic_ols | -4.06 | 6.58 | 6.17σ | 25,445 | OLS with full systematic model (cosD + cos2D + sin_m + cos_m + sin_y + cos_y + const) | SYSTEMATIC-CONTROLLED SENSITIVITY - upper bound; most leverage-sensitive |
| ar1_gls_full_model | -4.45 | 9.87 | 4.51σ | 25,445 | Full-model AR(1) GLS with Cochrane-Orcutt on full design matrix + cluster-robust SE (by station) | ROBUSTNESS CHECK - accounts for temporal autocorrelation with systematic controls |
| ar1_gls_cosd_only | -3.26 | 9.48 | 3.44σ | 25,445 | cosD-only AR(1) GLS with Cochrane-Orcutt + cluster-robust SE (by station) | COMPARISON - cosD-only; superseded by full-systematic model |
| full_sample_ols | -3.18 | 6.05 | 5.25σ | 25,445 | cosD-only OLS with 6σ-equivalent (MAD-based) outlier cleaning (step_003) | SECONDARY - cosD-only baseline |
| bayesian_mcmc | -2.87 | 6.63 | 4.32σ | 25,445 | Ensemble MCMC (32 walkers, 5000 steps) | SECONDARY - consistent with primary |
| leverage_excised_ols | -3.31 | 5.85 | 5.65σ | 25,176 | cosD-only OLS with Cook's Distance excision (threshold: 4/n) | DIAGNOSTIC - cosD-only leverage check |

### Robust Estimands

| Estimator | η (×10⁻⁴) | Error (×10⁻⁵) | SNR | Method | Status |
|-----------|-----------|----------------|-----|--------|--------|
| theil_sen | -2.94 | N/A | N/A | Median of pairwise slopes (cosD-only) | NONPARAMETRIC LOWER ENVELOPE - cosD-only; not directly comparable to full-systematic |
| precision_weighted_cosd_only | -3.21 | 5.12 | 6.27σ | cosD-only WLS with 1/σ² station weights | CROSS-STATION VALIDATION - cosD-only baseline |

### Table A: Station-level regression estimates

| Station | N | η (×10⁻⁴) | ση (×10⁻⁴) | η/ση |
|---------|---|-----------|------------|------|
| APO | 2,595 | -2.40 | 0.87 | 2.77 |
| Grasse | 18,742 | -3.78 | 0.66 | 5.73 |
| Haleakala | 666 | 25.71 | 11.65 | 2.21 |
| Matera | 345 | 3.85 | 8.30 | 0.46 |
| McDonald2 | 3,097 | -0.88 | 3.06 | 0.29 |

### Table B: Station-level phase-correlation diagnostics

| Station | r(R, cos D) | p-value | Phase coverage | RMS (cm) |
|---------|-------------|---------|----------------|----------|
| APO | -0.0543 | 5.68e-03 | Good | 3.16 |
| Grasse | -0.0418 | 9.96e-09 | Good | 5.90 |
| Haleakala | 0.0854 | 2.76e-02 | Biased: mean cos(D)=-0.347 | 10.62 |
| Matera | 0.0250 | 6.43e-01 | Biased: mean cos(D)=-0.025 | 5.89 |
| McDonald2 | -0.0052 | 7.74e-01 | Biased: mean cos(D)=-0.315 | 8.07 |