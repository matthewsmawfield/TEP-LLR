#!/usr/bin/env python3
"""
TEP-LLR: Full Canonical Analysis Pipeline
==========================================

Executes the complete 40-step LLR analysis pipeline in sequence.
Every step writes a JSON output to results/outputs/ and a detailed
log to logs/.  Steps are fail-fast: execution halts on the first
failure so that downstream steps do not consume stale data.

Pipeline Architecture
---------------------
Group A — Core Detection (steps 000–003)
    Ingest, preprocess, and perform the primary statistical detection
    of the cos(elongation) Nordtvedt signal.

Group B — Extended Validation & Robustness (steps 004–022)
    Multi-ephemeris comparison, systematic error characterisation,
    Bayesian inference, leverage diagnostics, temporal amplitude
    evolution, IPW station-balance validation, and 20 independent
    detection methods.

Group C — Physical Signal Probes (steps 023–028)
    Ephemeris-absorption masking simulation, environmental
    (heliocentric) amplitude modulation, solar-cycle correlation,
    thermal array deformation falsification, leverage temporal
    clustering, and TEP core-density suppression simulation.

Group D — Defensibility (steps 029–032)
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
        # 001b — Preprocess DE430 (JPL) residuals for cross-ephemeris comparison.
        "step_001b_de430_preprocessing.py",
        # 001c — Preprocess DE440 (JPL, 2021) residuals.
        "step_001c_de440_ephemeris.py",
        # 001d — Preprocess INPOP21a residuals.
        "step_001d_inpop21a_ephemeris.py",
        # 002 — Primary statistical analysis: Pearson r, OLS regression,
        #        differential (new moon vs full moon) analysis.  Produces the
        #        headline 7.9σ result.
        "step_002_statistical_analysis.py",
        # 003 — Advanced multi-method detection: bootstrap, permutation,
        #        Theil-Sen, leverage, outlier detection, station-by-station,
        #        temporal bins, phase bins, cross-validation, holdout test,
        #        Bayesian MCMC, Lomb-Scargle, Grand Phase Fold (20 methods).
        "step_003_detection_analysis_advanced.py",
        # -------------------------------------------------------------------
        # Group B: Extended Validation & Robustness
        # -------------------------------------------------------------------
        # 004 — Temporal drift analysis: linear and quadratic secular trends
        #        in eta over the 35-year baseline.
        "step_004_temporal_drift_analysis.py",
        # 005 — Multi-ephemeris comparison: compute eta from INPOP19a, DE430,
        #        DE440, INPOP21a in a unified framework.  Establishes cross-
        #        ephemeris consistency (both INPOP19a and cleaned DE430 are
        #        significant at >7σ).
        "step_005_multi_ephemeris_comparison.py",
        # 055 — Meta-analysis of ephemeris results: Bayesian combination of
        #        INPOP19a and DE430 with baseline-weighted uncertainties and
        #        systematic error quantification. Provides methodological
        #        strengthening through sign-consistency validation.
        "step_055_meta_analysis.py",
        # 006 — Construct systematic error budget (ephemeris, atmosphere,
        #        timing, tides).  Verify none correlate with cos(D).
        "step_006_systematic_error_analysis.py",
        # 007 — Ingest any additional station datasets not in the primary file.
        "step_007_additional_station_ingestion.py",
        # 008 — Apply systematic error corrections and produce corrected CSV.
        "step_008_systematic_error_correction.py",
        # 009 — Re-run primary analysis on corrected residuals; compare with
        #        uncorrected to confirm corrections don't suppress the signal.
        "step_009_corrected_data_analysis.py",
        # 010 — Ephemeris-independent analysis: test whether the signal
        #        survives when no ephemeris model parameters are assumed.
        "step_010_ephemeris_independent_analysis.py",
        # 011 — Systematic control analysis: partial correlations controlling
        #        for temporal trends, seasonal cycles, and station drifts.
        #        Signal persists at r=-0.034 (5.7σ) after all controls.
        "step_011_systematic_control_analysis.py",
        # 012 — Noise injection and signal recovery: signal survives 2×RMS
        #        noise addition; pipeline recovers injected signals at 100%.
        "step_012_noise_signal_injection.py",
        # 013 — Subsample robustness: five-category validation including
        #        station jackknife (all 4 powered leave-one-out samples
        #        significant), weight sensitivity, and IPW regression.
        "step_013_subsample_robustness.py",
        # 014 — Station decomposition: decompose the combined eta into station-
        #        specific contributions weighted by data fraction.
        "step_014_station_decomposition.py",
        # 014b — Inter-Station Consistency Analysis.
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
        #        Savage-Dickey Bayes Factor B = 1.8×10¹¹ (Decisive evidence).
        "step_016_bayesian_analysis.py",
        # 017 — TEP theoretical prediction: compute expected eta from the
        #        TEP suppression model and compare with measured value.
        "step_017_tep_prediction.py",
        # 018 — Leverage diagnostics: Cook's distance, hat-matrix, DFFITS.
        #        Identifies 5,016 high-leverage points (19.1%) whose presence
        #        accounts for the OLS/Theil-Sen factor-of-2 difference.
        "step_018_leverage_diagnostics.py",
        # 019 — Station quality metrics: per-station RMS, phase coverage,
        #        temporal baseline, detection power assessment.
        "step_019_station_quality.py",
        # 020 — Systematic Monte Carlo: 10,000 MC trials of plausible
        #        systematic scenarios; none reproduce the observed signal.
        "step_020_systematic_monte_carlo.py",
        # 021 — Temporal amplitude evolution: sliding-window eta(t) with
        #        trend decomposition.  Large chi²/dof reflects noise-level
        #        variation across hardware eras, not genuine signal variation.
        "step_021_temporal_amplitude.py",
        # 022 — IPW (Inverse-Probability Weighted) regression validation:
        #        equal per-station total weight.  Low IPW SNR is explained by
        #        phase-truncated McDonald2 (steps 031 diagnoses this fully).
        "step_022_ipw_validation.py",
        # -------------------------------------------------------------------
        # Group C: Physical Signal Probes
        # -------------------------------------------------------------------
        # 023 — Ephemeris absorption masking simulation: demonstrates that a
        #        dynamically oscillating eta is masked to near-zero by
        #        standard constant-parameter least-squares fitting, but
        #        survives intact in the post-fit residuals.
        "step_023_absorption_simulation.py",
        # 024 — Environmental modulation: test whether eta scales with
        #        heliocentric distance (perihelion vs aphelion split).
        #        TEP predicts weaker suppression at aphelion.
        "step_024_environmental_modulation.py",
        # 025 — Solar cycle correlation: test whether the sliding-window
        #        eta correlates with 11-year solar activity proxy.
        #        Permutation p=0.0007 confirms the correlation is non-random.
        "step_025_solar_cycle_correlation.py",
        # 026 — Thermal array deformation falsification: worst-case Apollo
        #        retroreflector thermal expansion is ~1 mm (11% of signal).
        "step_026_thermal_array_modeling.py",
        # 027 — Leverage temporal clustering: high-leverage points cluster
        #        in 1984–1989 (64.9% of Cook-flagged points from 8.4% of data),
        #        confirming early-era Nd:glass systematics drive OLS inflation.
        "step_027_leverage_temporal_clustering.py",
        # 028 — TEP core-density suppression simulation: derives theoretical
        #        eta from Earth's inner-core density gradient mapped against
        #        the Universal Critical Density, probing volumetric Temporal
        #        Topology flattening.
        "step_028_tep_core_density_simulation.py",
        # -------------------------------------------------------------------
        # Group D: Defensibility
        # -------------------------------------------------------------------
        # 029 — Day/night thermal bias null test: compute solar altitude for
        #        all 26,207 observations; partial regression shows the
        #        daytime-ranging bias explains <2.5% of the cos(D) amplitude.
        #        Result: daytime thermal hypothesis is rejected by the data.
        "step_029_day_night_thermal_bias.py",
        # 030 — True geometric elongation comparison: replace mean-phase
        #        cos(D_mean) with ephemeris-derived cos(D_true); partial
        #        regression shows D_mean dominates and D_true carries no
        #        independent signal.
        #        Result: mean-phase approximation is not a systematic bias.
        "step_030_geometric_elongation.py",
        # 031 — Station power analysis and Grasse-dominance assessment.
        #        Five-component analysis:
        #          (1) Per-station expected vs observed SNR
        #          (2) Phase-coverage quality (McDonald2 identified as biased:
        #              mean cos D = -0.326, only 3% of obs near full moon)
        #          (3) Monte Carlo IPW simulation
        #          (4) Precision-weighted (1/σ²) regression: eta=-4.94e-4 at 4.59σ
        #          (5) Grasse chronological split: both halves detect negative eta
        #        Status: CONSISTENT.
        "step_031_station_power_analysis.py",
        # 032 — Hardware epoch consistency analysis.
        #        Partitions data into 5 verified hardware epochs
        #        (Grasse Nd:glass/PMT 1984-93, Nd:YAG/SPAD 1994-2008,
        #        Nd:YAG/C-SPAD 2009-19; APO early 2000-09, mature 2010-19).
        #        All 5 epochs show negative eta; amplitude-RMS correlation
        #        r=0.988 (p=0.002) confirms scatter is noise-driven.
        #        Result: temporal non-stationarity is an instrumentation
        #        precision effect, not a sign-reversing systematic.
        "step_032_hardware_epoch_analysis.py",
        # 033 — Lomb-Scargle Orbital Dynamics Mapping.
        #        High-resolution periodogram mapping spectral residuals to
        #        formal Delaunay combinations, replacing ad-hoc approximations.
        "step_033_lomb_scargle_orbital_dynamics.py",
        # 034 — Ephemeris Orthogonality Proof.
        #        Demonstrates that a dynamically scaling TEP parameter
        #        structurally generates D +/- l' sideband harmonics, proving
        #        it is mathematically orthogonal to static Keplerian ephemeris bases.
        "step_034_ephemeris_orthogonality_proof.py",
        # 035 — Quantitative η Prediction from TEP First Principles.
        #        Derives predicted Nordtvedt parameter from TEP geometric
        #        suppression formalism and compares against measured η.
        "step_035_quantitative_eta_prediction.py",
        # 036 — Static vs Dynamic Signal Absorption Test.
        #        Validates that standard ephemeris solvers correctly recover
        #        static η but fail to absorb dynamic TEP sideband variance.
        "step_036_static_dynamic_absorption.py",
        # 037 — Historical Comparison Analysis.
        #        Quantitatively compares TEP detection to Müller & Nordtvedt (1998)
        #        unexplained ~1 cm synodic residual signal.
        "step_037_historical_comparison.py",
        # 038 — Full-Moon Deficit Analysis.
        #        Tests Murphy/Sabhlok full-moon deficit vs TEP scalar-field predictions.
        "step_038_full_moon_deficit_analysis.py",
        # 039 — Lunar Recession Analysis.
        #        Tests lunar orbit recession anomaly against TEP predictions.
        "step_039_lunar_recession_analysis.py",
        # 040 — Tidal Resonance Analysis.
        #        Tests North Atlantic tidal resonance explanation vs TEP alternative.
        "step_040_tidal_resonance_analysis.py",
        # 041 — Dust Model Sensitivity Analysis.
        #        Formal parameter sweep of Sabhlok et al. (2024) thermal/dust
        #        model showing dust estimate is underdetermined (20-80%) and
        #        thermal mechanism cannot explain observed 8.9 mm signal.
        "step_041_dust_sensitivity_analysis.py",
    ]

    run_pipeline("Full Canonical", steps, stop_on_failure=True)


if __name__ == "__main__":
    main()
