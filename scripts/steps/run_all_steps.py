#!/usr/bin/env python3
"""
TEP-LLR: Full Canonical Analysis Pipeline
==========================================

Executes the complete 44-step LLR analysis pipeline in sequence.
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
      029 — Day/Night thermal bias null test
      030 — True geometric elongation vs mean-phase comparison
      031 — Station power analysis and Grasse-dominance defense
             (includes phase-coverage analysis explaining McDonald2)
      032 — Hardware epoch consistency analysis (sign-consistency
             across Nd:glass / Nd:YAG / C-SPAD era partitions)

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
        # 000 — Download and parse INPOP19a LLR residual normal-point data
        #        from Paris Observatory.  Produces the master parquet file.
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
        # 014 — Inter-Station Consistency Analysis.
        #        Wald test, meta-analysis, and pairwise consistency across
        #        independent LLR stations (APO, Grasse, Matera, McDonald, Haleakala).
        "step_014_inter_station_consistency.py",
        # 015 — Null tests against control datasets (shuffled phases, DE430
        #        raw, non-TEP frequency bands) to calibrate false-alarm rates.
        "step_015_null_tests.py",
        # 016 — Bayesian MCMC analysis (emcee, 32 walkers × 10,000 steps).
        #        Initialised from fresh OLS estimate to ensure Gelman-Rubin
        #        convergence diagnostics are a genuine test of mixing
        #        (not biased by a hardcoded starting position).
        #        Savage-Dickey Bayes Factor B ≈ 3.8×10² (Strong evidence);
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
        #        TEP signal to assess out-of-sample predictive stability.
        "step_051_cross_validation_analysis.py",
        # 052 — Station Distribution Analysis.
        #        Quantifies the impact of station geographic distribution
        #        and observational sampling on detection robustness.
        "step_052_station_distribution_analysis.py",
        # 053 — Clean Subset Analysis.
        #        Repeats primary analysis on the highest-quality data subset
        #        to verify signal is not driven by low-quality observations.
        "step_053_clean_subset_analysis.py",
        # 054 — Toy Orbital TEP Perturbation.
        #        Analytic toy-model simulation of orbital perturbations
        #        induced by a temporally-varying Nordtvedt parameter.
        "step_054_toy_orbital_tep_perturbation.py",
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
        # 043 — Temporal Bin Variation Analysis.
        #        Quantitative analysis of temporal bin variation to address χ²/dof ≈ 33
        #        concern. Assesses whether temporal variation exceeds expected noise.
        "step_043_temporal_bin_variation_analysis.py",
        # 046 — Station-Balanced TEP Analysis.
        #        Directly addresses Grasse-dominance concern by creating balanced
        #        subsamples (equal-N per station, Grasse-capped). Tests whether
        #        signal persists when station contributions are equalised.
        "step_046_station_balanced_tep.py",
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
    ]

    run_pipeline("Full Canonical", steps, stop_on_failure=True)


if __name__ == "__main__":
    main()
