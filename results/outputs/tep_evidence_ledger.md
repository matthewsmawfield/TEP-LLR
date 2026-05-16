# TEP Evidence Ledger

Strong residual-channel evidence for TEP with bounded station and source-level-refit risks.

## Positive Evidence

### Primary synodic estimand (precision-weighted WLS)
- Source: `step_040_unified_results_table.json / step_050_corrected_tep_analysis.json`
- Interpretation: Headline Nordtvedt estimate from precision-weighted full-systematic WLS on all cleaned shots (no row deletion).
- eta: -3.909e-04
- eta_error: 5.628e-05
- snr: 6.94sigma
- snr_cluster: 6.78sigma
- n_obs: 25445

### Uncertainty calibration (prediction intervals and η brackets)
- Source: `step_064_prediction_coverage.json`
- Interpretation: Prediction intervals are conservative (observed coverage exceeds nominal at 68–95%). Published σ and/or the pooled regression variance scale are larger than residual scatter after the full-systematic fit (χ²_red < 1). Headline formal σ-based significance is therefore not inflated by underestimated errors; station-block bootstrap and LOSO conformal intervals provide σ-free significance brackets alongside WLS and cluster-robust paths.
- wls_68pct_observed: 8.591e-01
- wls_95pct_observed: 9.784e-01
- chi2_reduced_wls: 4.768e-01
- headline_snr_wls: 6.94sigma
- sigma_calibration_scale: 5.141e-01
- loso_conformal_95_excludes_zero: True

### Leverage diagnostic and cross-station support
- Source: `step_040_unified_results_table.json / step_050_corrected_tep_analysis.json`
- Interpretation: Cook's-Distance excised OLS and common-eta mixed model remain negative and significant; sign stable under leverage removal.
- cooks_excised_eta: -3.874e-04
- cooks_excised_snr: 7.82sigma
- cooks_excised_n_obs: 23837
- common_eta: -4.321e-04
- common_eta_snr: 6.40sigma
- station_deviation_p: 3.151e-01

### Known-systematic-only simulations
- Source: `step_061_systematic_sensitivity_analysis.json`
- Interpretation: No modeled known systematic at observed amplitude reproduces the full-systematic eta.
- required_amplitude_cm: 5.280e-01
- minimum_required_to_known_ratio: 1.486e+00
- minimum_ratio_systematic: ephemeris
- max_p_exceed_observed: 0.000e+00
- n_mc_per_systematic: 2000

### Adversarial nuisance + blind hold-out (Step 061)
- Source: `step_061_systematic_sensitivity_analysis.json`
- Interpretation: Data-driven PCA does not absorb cos(D); blind year hold-out recovers significant negative eta with nuisances trained only on non-held-out years.
- pca_joint_eta: -2.927e-04
- pca_joint_snr: 2.95sigma
- gp_absorption_fraction: 8.317e-01
- blind_holdout_eta: -3.658e-04
- blind_holdout_snr: 14.27sigma
- interaction_cells_fitted: 26
- interaction_cells_negative_eta: 13

### Linearized post-fit extraction
- Source: `step_056_dynamical_integrator_eta_refit.json`
- Interpretation: Published residual archives recover a consistent negative eta under the same full-systematic nuisance design.
- inpop_eta: -4.061e-04
- inpop_snr: 6.17sigma
- de430_eta: -5.977e-04
- de430_snr: 5.04sigma

### Directional residual structure
- Source: `step_055_cmb_rigorous_falsification.json`
- Interpretation: Corroborative fixed-sky directional anatomy on the residual channel, not a replacement eta estimator. Synodic phase coupling is decisively rejected; uniform-axis uniqueness remains marginal. Dual-axis fit: r(cos theta_CMB, cos theta_gal)=0.984; Planck term absorbed when both axes are included.
- diagnostics_passed: 4
- diagnostics_total: 5
- sky_scramble_p: 9.540e-02
- correlation_matched_p: 2.643e-01
- phase_null_p_eff: 6.000e-04
- orthogonal_scramble_p_eff: 9.020e-02
- dual_axis_r_cmb_gal: 9.844e-01
- dual_eta_cmb_t: 6.820e-01
- dual_eta_gal_t: -2.626e+00

### Leave-one-station-out meta (Step 072)
- Source: `step_072_leave_one_station_out_meta.json`
- Interpretation: Four of five LOSO exclusions remain powered with the same negative sign; excluding Grasse is underpowered. Inverse-variance meta on powered exclusions: 12.77sigma, I2=0%.
- eta_meta: -4.207e-04
- eta_meta_err: 3.294e-05
- snr_meta: 12.77sigma
- cochrans_Q: 5.040e-02
- I2_percent: 0.000e+00
- n_powered_loso: 4

## Bounded Risks

### Station leverage (LOSO + Grasse-conditioned)
- Source: `step_029_station_power_analysis.json / step_059_grasse_systematic_sufficiency.json / step_072_leave_one_station_out_meta.json / step_071_stratified_equal_n.json`
- Status: ['PASS', 'PASS', 'PASS', 'PASS']
- Interpretation: Grasse leverage is material: excluding Grasse is underpowered, but four other LOSO exclusions remain powered with the same negative sign (meta-analysis 12.77sigma). Non-Grasse and Grasse-conditioned estimands retain negative eta; equal-N balance reduces SNR without sign reversal.

### Source-level absorption/refit
- Source: `step_065_high_dimensional_absorption_test.json`
- Status: PASS
- Interpretation: High-dimensional residual-basis stress test supports sideband survival but preserves the need for source-level INPOP/DE430 eta-free refits.
