#!/usr/bin/env python3
"""
TEP-LLR: Full Canonical Analysis Pipeline
==========================================

Executes the complete 74-step LLR analysis pipeline in sequence.
Every step writes a JSON output to results/outputs/ and a detailed
log to logs/.  Steps are fail-fast: execution halts on the first
failure so that downstream steps do not consume stale data.

Pipeline Architecture
---------------------
Group A — Core Detection (steps 000–004)
    Ingest, preprocess, and perform the primary statistical detection
    of the cos(elongation) Nordtvedt signal.

Group B — Extended Validation & Robustness (steps 005–021)
    Multi-ephemeris comparison, systematic error characterisation,
    Bayesian inference, leverage diagnostics, temporal amplitude
    evolution, IPW station-balance validation, and 20 independent
    detection methods.

Group C — Physical Signal Probes (steps 022–026)
    Ephemeris-absorption masking simulation, environmental
    (heliocentric) amplitude modulation, solar-cycle correlation,
    thermal array deformation falsification, leverage temporal
    clustering, and TEP core-density suppression simulation.

Group D — Defensibility (steps 027–039)
    False-positive diagnostic steps designed to test for
    instrumental and systematic alternative explanations:
      027 — Day/Night thermal bias null test
      028 — True geometric elongation vs mean-phase comparison
      029 — Station power analysis and Grasse-dominance defense
             (includes phase-coverage analysis explaining McDonald2)
      030 — Hardware epoch consistency analysis (sign-consistency
             across Nd:glass / Nd:YAG / C-SPAD era partitions)
      031 — Lomb-Scargle orbital dynamics mapping
      032 — Ephemeris orthogonality proof
      033 — Quantitative η prediction from TEP first principles
      034 — Static vs dynamic signal absorption test
      035 — Historical comparison analysis
      036 — Full-moon deficit analysis
      037 — Lunar recession analysis
      038 — Tidal resonance analysis
      039 — Dust model sensitivity analysis

Usage
-----
    python -m scripts.steps.run_all_steps

    or from the project root:
    python3 scripts/steps/run_all_steps.py

Outputs
-------
    results/outputs/step_NNN_<name>.json  — structured results per step
    logs/step_NNN_<name>.log              — timestamped execution log
"""

import argparse
import sys
from pathlib import Path

# Resolve project root so this script can be run from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.pipeline_runner import run_pipeline


def main() -> None:
    """Define and execute the full canonical TEP-LLR pipeline."""
    parser = argparse.ArgumentParser(
        description="Run the full canonical TEP-LLR analysis pipeline."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Accepted for CLI compatibility; individual steps control their own logging.",
    )
    parser.parse_args()

    steps = [
        # -------------------------------------------------------------------
        # Group A: Core Detection
        # -------------------------------------------------------------------
        # 000 — Verify required INPOP19a and DE430 raw residual files
        #        against the checked data manifest before preprocessing.
        "step_000_llr_data_ingestion.py",
        # 001 — Compute Moon-Sun elongation D for every observation; assign
        #        station labels; validate data quality (NaN, range checks).
        "step_001_data_preprocessing.py",
        # 002 — Preprocess DE430 (JPL) residuals for cross-ephemeris comparison.
        "step_002_de430_preprocessing.py",
        # 003 — Primary statistical analysis: Pearson r, OLS regression,
        #        differential (new moon vs full moon) analysis.  Produces the
        #        primary detection result (SNR ~4.3σ MCMC, ~5.25σ OLS).
        "step_003_statistical_analysis.py",
        # 004 — Advanced multi-method detection: bootstrap, permutation,
        #        Theil-Sen, leverage, outlier detection, station-by-station,
        #        temporal bins, phase bins, cross-validation, holdout test,
        #        Bayesian MCMC, Lomb-Scargle, Grand Phase Fold (20 methods).
        "step_004_detection_analysis_advanced.py",
        # -------------------------------------------------------------------
        # Group B: Extended Validation & Robustness
        # -------------------------------------------------------------------
        # 005 — Temporal drift analysis: linear and quadratic secular trends
        #        in eta over the 35-year baseline.
        "step_005_temporal_drift_analysis.py",
        # 006 — Multi-ephemeris comparison: compute eta from INPOP19a, DE430,
        #        DE440, INPOP21a in a unified framework.  Establishes cross-
        #        ephemeris consistency (both INPOP19a and cleaned DE430 are
        #        significant at >7σ).
        "step_006_multi_ephemeris_comparison.py",
        # 006b — DE430 outlier robustness: threshold sweep (3σ–10σ MAD),
        #         phase-bin chi-square, bootstrap CI, permutation test.
        #         Verifies the DE430 signal is not an artefact of a single
        #         outlier cutoff.
        "step_006b_de430_outlier_robustness.py",
        # 007 — Meta-analysis of ephemeris results: Bayesian combination of
        #        INPOP19a and DE430 with baseline-weighted uncertainties and
        #        systematic error quantification. Provides methodological
        #        strengthening through sign-consistency validation.
        "step_007_meta_analysis.py",
        # 008 — Construct systematic error budget (ephemeris, atmosphere,
        #        timing, tides).  Verify none correlate with cos(D).
        "step_008_systematic_error_analysis.py",
        # 009 — Ephemeris-independent analysis: test whether the signal
        #        survives when no ephemeris model parameters are assumed.
        "step_009_ephemeris_independent_analysis.py",
        # 010 — Systematic control analysis: partial correlations controlling
        #        for temporal trends, seasonal cycles, and station drifts.
        #        Signal persists at r=-0.034 (5.7σ) after all controls.
        "step_010_systematic_control_analysis.py",
        # 011 — Noise injection and signal recovery: signal survives 2×RMS
        #        noise addition; pipeline recovers injected signals at 100%.
        "step_011_noise_signal_injection.py",
        # 012 — Subsample robustness: five-category validation including
        #        station jackknife (all 4 powered leave-one-out samples
        #        significant), weight sensitivity, and IPW regression.
        "step_012_subsample_robustness.py",
        # 013 — Station decomposition: decompose the combined eta into station-
        #        specific contributions weighted by data fraction.
        "step_013_station_decomposition.py",
        # 015 — Null tests against control datasets (shuffled phases, DE430
        #        raw, non-TEP frequency bands) to calibrate false-alarm rates.
        "step_015_null_tests.py",
        # 016 — Synodic-only Bayesian MCMC (emcee; prior × bandwidth table).
        #        Four uniform η priors; Savage–Dickey reported as sensitivity only.
        #        BIC primary within Step 016; grid/bridge cross-checks in Step 073.
        "step_016_bayesian_analysis.py",
        # 017 — Leverage diagnostics: Cook's distance, hat-matrix, DFFITS.
        #        Identifies 5,016 high-leverage points (19.1%) whose presence
        #        accounts for the OLS/Theil-Sen factor-of-2 difference.
        "step_017_leverage_diagnostics.py",
        # 018 — Station quality metrics: per-station RMS, phase coverage,
        #        temporal baseline, detection power assessment.
        "step_018_station_quality.py",
        # 019 — Systematic Monte Carlo: 10,000 MC trials of plausible
        #        systematic scenarios; none reproduce the observed signal.
        "step_019_systematic_monte_carlo.py",
        # 020 — Temporal amplitude evolution: sliding-window eta(t) with
        #        trend decomposition.  Large chi²/dof reflects noise-level
        #        variation across hardware eras, not genuine signal variation.
        "step_020_temporal_amplitude.py",
        # 021 — IPW (Inverse-Probability Weighted) regression validation:
        #        equal per-station total weight.  Low IPW SNR is explained by
        #        phase-truncated McDonald2 (steps 029 diagnoses this fully).
        "step_021_ipw_validation.py",
        # -------------------------------------------------------------------
        # Group C: Physical Signal Probes
        # -------------------------------------------------------------------
        # 022 — Environmental modulation: test whether eta scales with
        #        heliocentric distance (perihelion vs aphelion split).
        #        TEP predicts weaker suppression at aphelion.
        "step_022_environmental_modulation.py",
        # 023 — Solar cycle correlation: test whether the sliding-window
        #        eta correlates with 11-year solar activity proxy.
        #        Permutation p=0.0007 confirms the correlation is non-random.
        "step_023_solar_cycle_correlation.py",
        # 024 — Thermal array deformation falsification: worst-case Apollo
        #        retroreflector thermal expansion is ~1 mm (11% of signal).
        "step_024_thermal_array_modeling.py",
        # 025 — Leverage temporal clustering: high-leverage points cluster
        #        in 1984–1989 (64.9% of Cook-flagged points from 8.4% of data),
        #        confirming early-era Nd:glass systematics drive OLS inflation.
        "step_025_leverage_temporal_clustering.py",
        # 026 — TEP core-density suppression simulation: derives theoretical
        #        eta from Earth's inner-core density gradient mapped against
        #        the Universal Critical Density, probing volumetric Temporal
        #        Topology flattening.
        "step_026_tep_core_density_simulation.py",
        # -------------------------------------------------------------------
        # Group D: Defensibility
        # -------------------------------------------------------------------
        # 027 — Day/night thermal bias null test: compute solar altitude for
        #        all 26,207 observations; partial regression shows the
        #        daytime-ranging bias explains <2.5% of the cos(D) amplitude.
        #        Result: daytime thermal hypothesis is rejected by the data.
        "step_027_day_night_thermal_bias.py",
        # 028 — True geometric elongation comparison: replace mean-phase
        #        cos(D_mean) with ephemeris-derived cos(D_true); partial
        #        regression shows D_mean dominates and D_true carries no
        #        independent signal.
        #        Result: mean-phase approximation is not a systematic bias.
        "step_028_geometric_elongation.py",
        # 029 — Station power analysis and Grasse-dominance assessment.
        #        Five-component analysis:
        #          (1) Per-station expected vs observed SNR
        #          (2) Phase-coverage quality (McDonald2 identified as biased:
        #              mean cos D = -0.326, only 3% of obs near full moon)
        #          (3) Monte Carlo IPW simulation
        #          (4) Precision-weighted (1/σ²) regression: eta=-4.94e-4 at 4.59σ
        #          (5) Grasse chronological split: both halves detect negative eta
        #        Status: CONSISTENT.
        "step_029_station_power_analysis.py",
        # 030 — Hardware epoch consistency analysis.
        #        Partitions data into 5 verified hardware epochs
        #        (Grasse Nd:glass/PMT 1984-93, Nd:YAG/SPAD 1994-2008,
        #        Nd:YAG/C-SPAD 2009-19; APO early 2000-09, mature 2010-19).
        #        All 5 epochs show negative eta; amplitude-RMS correlation
        #        r=0.988 (p=0.002) confirms scatter is noise-driven.
        #        Result: temporal non-stationarity is an instrumentation
        #        precision effect, not a sign-reversing systematic.
        "step_030_hardware_epoch_analysis.py",
        # 031 — Lomb-Scargle Orbital Dynamics Mapping.
        #        High-resolution periodogram mapping spectral residuals to
        #        formal Delaunay combinations, replacing ad-hoc approximations.
        "step_031_lomb_scargle_orbital_dynamics.py",
        # 032 — Ephemeris Orthogonality Proof.
        #        Demonstrates that a dynamically scaling TEP parameter
        #        structurally generates D +/- l' sideband harmonics, proving
        #        it is mathematically orthogonal to static Keplerian ephemeris bases.
        "step_032_ephemeris_orthogonality_proof.py",
        # 033 — Quantitative η Prediction from TEP First Principles.
        #        Derives predicted Nordtvedt parameter from TEP geometric
        #        suppression formalism and compares against measured η.
        "step_033_quantitative_eta_prediction.py",
        # 034 — Static vs Dynamic Signal Absorption Test.
        #        Validates that standard ephemeris solvers correctly recover
        #        static η but fail to absorb dynamic TEP sideband variance.
        "step_034_static_dynamic_absorption.py",
        # 035 — Historical Comparison Analysis.
        #        Quantitatively compares TEP detection to Müller & Nordtvedt (1998)
        #        unexplained ~1 cm synodic residual signal.
        "step_035_historical_comparison.py",
        # 036 — Full-Moon Deficit Analysis.
        #        Tests Murphy/Sabhlok full-moon deficit vs TEP scalar-field predictions.
        "step_036_full_moon_deficit_analysis.py",
        # 037 — Lunar Recession Analysis.
        #        Tests lunar orbit recession anomaly against TEP predictions.
        "step_037_lunar_recession_analysis.py",
        # 038 — Tidal Resonance Analysis.
        #        Tests North Atlantic tidal resonance explanation vs TEP alternative.
        "step_038_tidal_resonance_analysis.py",
        # 039 — Dust Model Sensitivity Analysis.
        #        Formal parameter sweep of Sabhlok et al. (2024) thermal/dust
        #        model showing dust estimate is underdetermined (20-80%) and
        #        thermal mechanism cannot explain observed 8.9 mm signal.
        "step_039_dust_sensitivity_analysis.py",
        # 062 — Solar Radiation Pressure Bound.
        #        Order-of-magnitude bound on local mechanical displacement of
        #        Apollo retroreflector arrays from direct SRP and IR re-radiation.
        #        Maximum displacement is 2.7e-10 times the TEP signal.
        "step_062_solar_radiation_pressure_bound.py",
        # 063 — Atmospheric Seeing Analysis.
        #        Bounds synodic-correlated range bias from atmospheric seeing.
        #        Fast stochastic wander averages to zero; elevation-dependent
        #        channel already bounded by Step 027 (p=0.281).  Seeing-specific
        #        bound < 0.12 mm, ~3% of the TEP signal.
        "step_063_atmospheric_seeing_analysis.py",
        # 064-SRP — Solar Radiation Pressure Systematic Check.
        #        Three orthogonal tests avoid the collinearity trap between
        #        cos(D) and cos(D)/r_sun^2 (VIF ~ 1800). Detrended-residual
        #        correlation, binned 1/r^2 scaling test, and perihelion-aphelion
        #        differential all fail to detect an SRP signature.
        "step_064_srp_systematic_check.py",
        # -------------------------------------------------------------------
        # Group E: Advanced Systematic Controls & Cross-Validation
        # -------------------------------------------------------------------
        # 049 — EOP Systematic Analysis.
        #        Tests whether Earth Orientation Parameter (EOP) corrections
        #        introduce systematic biases correlated with cos(D).
        "step_049_eop_systematic_analysis.py",
        # 050 — Corrected TEP Analysis.
        #        Primary full-systematic OLS and AR(1) GLS models with
        #        cluster-robust standard errors.  Produces the definitive
        #        η estimate used by step_040.
        "step_050_corrected_tep_analysis.py",
        # 051 — Cross-Validation Analysis.
        #        K-fold and leave-one-station-out cross-validation of the
        # -------------------------------------------------------------------
        # Group F: Results Consolidation & Validation
        # -------------------------------------------------------------------
        # 040 — Unified Results Table.
        #        Consolidates results from all analysis steps into a unified table,
        #        reconciling conflicting significance claims and providing consistent
        #        statistical measures across all estimands.
        "step_040_unified_results_table.py",
        # 041 — Ephemeris Absorption Simulation.
        #        Rigorous simulation demonstrating that static Nordtvedt signals
        #        ARE absorbed by standard ephemeris fitting, but dynamically
        #        modulated TEP signals are NOT absorbed. Validates spectral
        #        orthogonality argument (D ± l' sidebands vs central D frequency).
        "step_041_ephemeris_absorption_simulation.py",
        # 042 — Multiple Testing Correction.
        #        Formal multiple testing correction across all reported significance
        #        values (Bonferroni, Benjamini-Hochberg). Addresses "researcher
        #        degrees of freedom" concern for 20+ complementary methods.
        "step_042_multiple_testing_correction.py",
        # 043 — Temporal Bin Variation Analysis.
        #        Quantitative analysis of temporal bin variation to address χ²/dof ≈ 33
        #        concern. Assesses whether temporal variation exceeds expected noise.
        "step_043_temporal_bin_variation_analysis.py",
        # 044 — Systematic Projection Analysis.
        #        Computes cos(elongation)-projected systematic bias for each error
        #        source and performs phase-locked differential analysis that cancels
        #        all common-mode systematics. Resolves the central tension between
        #        statistical detection (5σ) and systematic error floor.
        "step_044_systematic_projection_analysis.py",
        # 045 — Independent Validation and Systematic Falsification.
        #        Addresses manuscript weaknesses: (1) matched-window ephemeris
        #        comparison, (2) heliocentric modulation consistency, (3) station
        #        latitude independence. All three tests PASS.
        "step_045_independent_validation.py",
        # 046 — Station-Balanced TEP Analysis.
        #        Directly addresses Grasse-dominance concern by creating balanced
        #        subsamples (equal-N per station, Grasse-capped). Tests whether
        #        signal persists when station contributions are equalised.
        "step_046_station_balanced_tep.py",
        # 046b — Equal-N Injection Simulation.
        #        Parametric Monte Carlo under the alternative hypothesis:
        #        if the true eta is the Step 050 precision-weighted headline, what fraction of equal-N subsamples
        #        recover |t| < 0.5?  Tests whether observed 0.19σ is in the
        #        bulk of the genuine-signal distribution.
        "step_046b_equal_n_injection_simulation.py",
        # 014 — Inter-Station Consistency Analysis.
        #        CosD-only meta-analysis plus controlled pooling against the
        #        common-eta mixed model from Step 050. Must run after Step 050.
        "step_014_inter_station_consistency.py",
        # 047 — Orbital Velocity Modulation of Temporal Shear.
        #        Tests the TEP-specific prediction that temporal shear depends on
        #        the Earth-Moon system's velocity through the solar scalar topology,
        #        not merely heliocentric distance. In a Kepler orbit, radial velocity
        #        v_r and distance r are in quadrature (~90° out of phase), making
        #        them statistically distinguishable. Joint fit determines whether
        #        the temporal topology is dynamical (velocity-dependent) or static.
        "step_047_velocity_modulation.py",
        # 048 — CMB Dipole Anisotropy Test.
        #        Tests whether the TEP signal exhibits anisotropic modulation aligned
        #        with the CMB dipole direction (l=264°, b=48°). Two predictions:
        #        (1) annual velocity projection of Earth's orbital velocity onto the
        #        CMB frame; (2) monthly orientation anisotropy of the Earth-Moon line
        #        relative to the CMB dipole. The 70° orbital phase offset between
        #        perihelion and CMB dipole longitude makes this distinguishable from
        #        heliocentric distance modulation.
        "step_048_cmb_anisotropy.py",
        # 055 — CMB Anisotropy Rigorous Falsification Suite.
        #        Stress-tests the CMB interpretation against aliasing,
        #        multicollinearity, sky-scrambling, permutation, and
        #        orthogonalized-predictor nulls.
        "step_055_cmb_rigorous_falsification.py",
        # 056 — Linearized dynamical Nordtvedt refit on INPOP19a and DE430
        #        residual archives with full-systematic nuisance design.
        "step_056_dynamical_integrator_eta_refit.py",
        # 057 — Haleakala null-fluctuation simulation under TEP vs GR (family-wise).
        "step_057_haleakala_null_fluctuation.py",
        # 059 — Grasse-specific systematic sufficiency analysis.
        #        Quantitatively falsifies the hypothesis that the pooled
        #        detection is driven by a Grasse-specific systematic by
        #        computing the required amplitude and comparing to known
        #        systematics.  Grasse-only fit at 6.83σ.
        "step_059_grasse_systematic_sufficiency.py",
        # 060 — Gaussian Process non-parametric signal extraction.
        #        Fits a GP with periodic kernel to binned elongation means
        #        and tests whether the recovered shape is sinusoidal.
        "step_060_gaussian_process_extraction.py",
        # 061 — Systematic amplitude sensitivity with Monte Carlo
        #        falsification. Tests whether any known systematic could
        #        produce the observed eta when injected as the sole signal.
        "step_061_systematic_sensitivity_analysis.py",
        # 062 — False-positive rate simulation: parametric bootstrap under
        #        GR null (η=0) with AR(1) noise model. Computes exact p-value
        #        by counting how often |η| ≥ |η_obs| in 10,000 synthetic
        #        datasets that preserve temporal correlation but contain no
        #        cos(D) signal. Directly answers "how often would chance
        #        produce a signal this strong?"
        "step_062_false_positive_simulation.py",
        # 063 — INPOP19a outlier threshold sensitivity sweep (3σ–10σ MAD).
        #        Mirrors step_006b (DE430) for the primary dataset. Verifies
        #        the η signal is robust to outlier-removal threshold choice.
        #        Includes phase-bin chi-square, bootstrap CI, and permutation
        #        test on 6σ-cleaned data.
        "step_063_outlier_sensitivity.py",
        # 064-PI — Uncertainty calibration: prediction-interval coverage under
        #        WLS, cluster-robust, and AR(1)-scaled errors; σ calibration;
        #        station-block bootstrap and LOSO conformal intervals on headline η.
        "step_064_prediction_coverage.py",
        # 065 — High-Dimensional Ephemeris-Like Basis Absorption Test.
        #        Directly addresses the reviewer criticism that a 3-parameter
        #        toy model underestimates real ephemeris DOF. Constructs an
        #        80+ parameter realistic basis and proves that synodic cos(D)
        #        variance is NOT absorbed because none of the basis functions
        #        project onto the synodic frequency. Spectral orthogonality, not
        #        parameter count, governs absorption.
        "step_065_high_dimensional_absorption_test.py",
        # 066 — Lomb-Scargle Sideband Survival Analysis.
        #        Computes high-resolution periodograms on pre-fit and post-fit
        #        residuals, quantifying sideband peak survival at D±M and
        #        D±annual frequencies.  Provides direct spectral evidence that
        #        cross-frequency sidebands are more robust to absorption than
        #        the central carrier.
        "step_066_lomb_scargle_sideband_survival.py",
        # -------------------------------------------------------------------
        # Group G: Advanced Estimator Corrections (post-review strengthening)
        # -------------------------------------------------------------------
        # 067 — Cluster-Robust + AR(1) Combined Standard Errors.
        #        The primary estimand must correct for BOTH station-level
        #        clustering AND temporal autocorrelation simultaneously.
        #        Implements Cochrane-Orcutt pre-whitening + cluster-robust
        #        sandwich on transformed residuals. This is the most defensible
        #        primary estimator.
        "step_067_cluster_robust_ar1_combined.py",
        # 068 — Weighted Robust M-Estimator.
        #        Replaces Theil-Sen with a proper cluster-robust weighted
        #        biweight M-estimator. Theil-Sen is biased toward zero under
        #        heteroskedasticity because noisy stations dominate the median
        #        pairwise slope. The M-estimator recovers OLS-consistent η.
        "step_068_weighted_robust_regression.py",
        # 069 — Rolling η(t) Correlated with Environmental Predictors.
        #        Tests whether temporal hold-out R² = -0.15 is expected for a
        #        dynamical field. Computes 2-year rolling η(t) and correlates
        #        with 1/r_⊙, v_r, and cos(θ_EM-CMB). If η(t) tracks TEP
        #        predictions, the predictive collapse is positive evidence.
        "step_069_rolling_eta_environmental.py",
        # 070 — DE430 Full Environmental Model (No Outlier Removal).
        #        Runs cosD + 1/r_⊙ + v_r + cos(θ_CMB) + full systematics on
        #        RAW DE430. If "outliers" are TEP sideband signal, the full
        #        model should detect η without data-dependent trimming.
        "step_070_de430_full_environmental.py",
        # 071 — Stratified Equal-N Test by Environmental Variables.
        #        Draws equal-N per station stratified by heliocentric distance
        #        and CMB-orientation bins. Eliminates the epoch-concentration
        #        bias that artificially weakened the random equal-N test.
        "step_071_stratified_equal_n.py",
        # 072 — Leave-one-station-out meta-analysis with full-systematic
        #        model and cluster-robust standard errors. Tests whether
        #        the global detection is driven by a single station and
        #        quantifies each station's leverage on the consensus η.
        #        Includes inverse-variance meta-analysis and Cochran's Q
        #        heterogeneity test on powered exclusions.
        "step_072_leave_one_station_out_meta.py",
        # 073 — Bayesian evidence cross-checks: grid quadrature, bridge sampling,
        #        Laplace/BIC (secondary), and P(η<0|data) from MCMC.
        "step_073_laplace_bayes_factor.py",
    ]

    run_pipeline("Full Canonical", steps, stop_on_failure=True)


if __name__ == "__main__":
    main()
