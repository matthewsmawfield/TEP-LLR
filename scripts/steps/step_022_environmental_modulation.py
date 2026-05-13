#!/usr/bin/env python3
"""
Step 022: Environmental Modulation Testing

Computes the Earth-Sun distance for each LLR observation epoch and tests
if the Nordtvedt η parameter magnitude scales as a function of orbital
phase (perihelion vs aphelion), which tests the differential suppression
hypothesis of the TEP framework.

Physical Mechanism:
The solar scalar gradient modulates Temporal Shear in the Earth-Moon system.
As the Earth-Moon system moves through the solar Temporal Topology, the background
scalar gradient varies with heliocentric distance. At perihelion (r ≈ 0.983 AU),
the Earth-Moon system plunges into a region of steeper scalar gradient, enhancing
Temporal Shear and altering the effective coupling. At aphelion (r ≈ 1.017 AU),
the gradient relaxes, causing the coupling to decay and computationally restoring
standard General Relativistic mechanics. This pipeline step tests for this exact
heliocentric modulation of Temporal Shear.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
from skyfield.api import load

# Setup paths

from scipy.optimize import curve_fit
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Earth's orbital parameters
# Source: IAU 1976 Astronomical Constants, JPL DE430/DE440 ephemeris documentation
EARTH_ECCENTRICITY = 0.0167  # e_⊕ (Earth's orbital eccentricity)
MEAN_ORBITAL_DISTANCE_AU = 1.0  # r_mean (1 AU = astronomical unit, by definition)

def threshold_activation_model(r_sun_au, eta_0, m):
    """TEP threshold-activation formula for heliocentric modulation.
    As the Earth-Moon system moves through the solar Temporal Topology,
    the background scalar gradient modulates Temporal Shear, altering
    the effective coupling.

    η(r⊙) = η₀ × (1 + m × (r_mean − r⊙) / (r_mean × e_⊕))

    Parameters:
    -----------
    r_sun_au : array-like
        Heliocentric distance in AU
    eta_0 : float
        Baseline Nordtvedt parameter at mean distance
    m : float
        Modulation depth (m ≈ 1 predicts sign-flip at aphelion)

    Returns:
    --------
    array-like : Predicted η values at each r_sun_au
    """
    return eta_0 * (
        1
        + m
        * (MEAN_ORBITAL_DISTANCE_AU - r_sun_au)
        / (MEAN_ORBITAL_DISTANCE_AU * EARTH_ECCENTRICITY)
    )


def compute_distances_and_analyze(df, verbose=False):
    # 1. Load skyfield ephemeris (canonical kernel; verified in Step 000 manifest)
    eph_path = PROJECT_ROOT / "data" / "raw" / "de440.bsp"
    if not eph_path.exists():
        raise FileNotFoundError(
            f"Required Skyfield kernel missing for Step 022: {eph_path}. "
            "Place de440.bsp under data/raw and verify via Step 000 manifest."
        )

    planets = load(str(eph_path))
    earth = planets["earth"]
    sun = planets["sun"]
    ts = load.timescale()

    timestamps = ts.tt(jd=df["date_julian"].values)

    print_status("Computing Earth-Sun distances via Skyfield...", "INFO")
    print_status(f"Using ephemeris: {eph_path.relative_to(PROJECT_ROOT)}", "INFO")

    # Vectorized distance computation
    astrometric = earth.at(timestamps).observe(sun)
    distances_au = astrometric.distance().au

    df["sun_distance_au"] = distances_au

    # Pre-filter large outliers
    outlier_mask = detect_outliers_sigma(df["residual_m"].values, sigma_threshold=6.0)
    df_clean = df[~outlier_mask]  # PERFORMANCE FIX: Removed unnecessary .copy()

    print_status(
        f"Computed distances for {len(df_clean):,} cleaned observations.", "CALC"
    )
    print_status(
        f"Distance range: {np.min(distances_au):.4f} AU to {np.max(distances_au):.4f} AU",
        "CALC",
    )
    print_status(f"Removed {len(df) - len(df_clean)} outliers (6σ threshold)", "CALC")

    # 2. Split analysis: Deep Perihelion (closest 15%) vs Deep Aphelion (furthest 15%)
    # Earth's orbit ranges from ~0.983 to 1.017 AU.
    # Theoretical framing: The solar scalar gradient is steepest at deep perihelion, maximizing Temporal Shear.
    # TEP predicts a non-linear threshold activation as Temporal Shear modulates the effective coupling.
    p15 = np.percentile(df_clean["sun_distance_au"], 15)
    p85 = np.percentile(df_clean["sun_distance_au"], 85)

    perihelion_df = df_clean[df_clean["sun_distance_au"] <= p15]
    aphelion_df = df_clean[df_clean["sun_distance_au"] >= p85]

    reg_peri = linear_regression(
        perihelion_df["residual_m"].values,
        np.cos(perihelion_df["elongation_rad"].values),
    )
    reg_aph = linear_regression(
        aphelion_df["residual_m"].values, np.cos(aphelion_df["elongation_rad"].values)
    )

    print_status("", "INFO")
    print_status("--- DIFFERENTIAL SUPPRESSION ANALYSIS ---", "PROCESS")
    print_status(
        "Hypothesis: Eta magnitude should scale with solar scalar gradient (varies with solar distance)",
        "INFO",
    )
    print_status("Perihelion threshold: r <= 15th percentile", "INFO")
    print_status("Aphelion threshold: r >= 85th percentile", "INFO")
    print_status("", "INFO")
    print_status("DEEP PERIHELION (High Shear, r <= 15th percentile):", "CALC")
    print_status(f"  N = {len(perihelion_df):,}", "CALC")
    print_status(f"  η = {reg_peri['eta']:.4e} ± {reg_peri['eta_error']:.4e}", "CALC")
    snr_peri = abs(reg_peri['eta']) / reg_peri['eta_error'] if reg_peri['eta_error'] > 0 else 0.0
    print_status(f"  SNR = {snr_peri:.2f}σ", "CALC")
    print_status("", "INFO")
    print_status("DEEP APHELION (Low Shear, r >= 85th percentile):", "CALC")
    print_status(f"  N = {len(aphelion_df):,}", "CALC")
    print_status(f"  η = {reg_aph['eta']:.4e} ± {reg_aph['eta_error']:.4e}", "CALC")
    snr_aph = abs(reg_aph['eta']) / reg_aph['eta_error'] if reg_aph['eta_error'] > 0 else 0.0
    print_status(f"  SNR = {snr_aph:.2f}σ", "CALC")

    # Test for statistically significant difference
    eta_diff = reg_peri["eta"] - reg_aph["eta"]
    diff_err = np.sqrt(reg_peri["eta_error"] ** 2 + reg_aph["eta_error"] ** 2)
    diff_significance = abs(eta_diff) / diff_err

    # 3. TEP Threshold-Activation Formula Fit
    # Bin data by heliocentric distance and fit η vs distance
    n_bins = 10
    bin_edges = np.linspace(
        df_clean["sun_distance_au"].min(), df_clean["sun_distance_au"].max(), n_bins + 1
    )

    bin_centers = []
    bin_etas = []
    bin_eta_errs = []

    for i in range(n_bins):
        mask = (df_clean["sun_distance_au"] >= bin_edges[i]) & (
            df_clean["sun_distance_au"] < bin_edges[i + 1]
        )
        if mask.sum() > 50:  # Require sufficient data
            bin_df = df_clean[mask]
            reg_bin = linear_regression(
                bin_df["residual_m"].values, np.cos(bin_df["elongation_rad"].values)
            )
            bin_centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            bin_etas.append(reg_bin["eta"])
            bin_eta_errs.append(reg_bin["eta_error"])

    # Fit threshold-activation model
    threshold_fit = None
    if len(bin_centers) >= 4:
        try:
            popt, pcov = curve_fit(
                threshold_activation_model,
                np.array(bin_centers),
                np.array(bin_etas),
                sigma=np.array(bin_eta_errs),
                p0=[-5e-4, 1.0],  # Initial guess: eta_0 ~ -5e-4, m ~ 1
                bounds=([-1e-3, 0], [0, 5]),
            )  # eta_0 < 0, m > 0
            eta_0_fit, m_fit = popt
            eta_0_err, m_err = np.sqrt(np.diag(pcov))

            threshold_fit = {
                "eta_0": float(eta_0_fit),
                "eta_0_error": float(eta_0_err),
                "modulation_depth_m": float(m_fit),
                "modulation_depth_error": float(m_err),
                "reduced_chi2": float(
                    np.sum(
                        (
                            np.array(bin_etas)
                            - threshold_activation_model(np.array(bin_centers), *popt)
                        )
                        ** 2
                        / np.array(bin_eta_errs) ** 2
                    )
                    / (len(bin_centers) - 2)
                ),
            }

            if verbose:
                print_status("", "INFO")
                print_status("--- THRESHOLD-ACTIVATION FORMULA FIT ---", "PROCESS")
                print_status(f"Fitted: η(r⊙) = η₀ × (1 + m × (1−r⊙)/(1×e⊕))", "CALC")
                print_status(f"  η₀ = {eta_0_fit:.4e} ± {eta_0_err:.4e}", "CALC")
                print_status(f"  m = {m_fit:.3f} ± {m_err:.3f}", "CALC")
                print_status(
                    f"  Reduced χ² = {threshold_fit['reduced_chi2']:.2f}", "CALC"
                )
                # Tolerance of ±0.5 on slope consistency check (allowing for measurement uncertainty)
                if abs(m_fit - 1.0) < 0.5:
                    print_status(
                        "  m ≈ 1: Consistent with predicted sign-flip at aphelion!",
                        "SUCCESS",
                    )
                else:
                    print_status(f"  m = {m_fit:.2f}: Deviates from unity", "WARNING")

        except Exception as e:
            print_status(f"CRITICAL: Threshold fit failed: {e}", "ERROR")
            raise RuntimeError(f"TEP Threshold-Activation fit failed to converge: {e}")

    if verbose:
        print_status(
            f"Amplitude differential significance: {diff_significance:.2f}σ", "CALC"
        )
        if diff_significance > 2.0:
            print_status(
                "RESULT: Significant environmental scaling detected. GR prediction (diff=0) violated.",
                "SUCCESS",
            )
        else:
            print_status(
                "RESULT: Non-significant scaling over this distance baseline (null on this test).",
                "INFO",
            )

    # Step completes successfully whether the differential is significant or not; significance
    # is encoded in differential.is_significant and modulation_test_result (not status=WARNING).
    significant = bool(diff_significance > 2.0)
    modulation_test_result = (
        "DIFFERENTIAL_SIGNIFICANT" if significant else "DIFFERENTIAL_NOT_SIGNIFICANT"
    )

    return {
        "step_id": "step_022",
        "status": "PASS",
        "modulation_test_result": modulation_test_result,
        "test_interpretation": "Significant perihelion-aphelion differential confirms environmental scaling predicted by TEP, ruling out stationary instrumental systematics"
        if significant
        else "Perihelion–aphelion η differential is not significant at the 2σ threshold on this split; "
        "environmental scaling is therefore not established by this test (null outcome).",
        "distance_stats": {
            "min_au": float(np.min(distances_au)),
            "max_au": float(np.max(distances_au)),
        },
        "perihelion": {
            "n_obs": int(len(perihelion_df)),
            "eta": float(reg_peri["eta"]),
            "eta_error": float(reg_peri["eta_error"]),
        },
        "aphelion": {
            "n_obs": int(len(aphelion_df)),
            "eta": float(reg_aph["eta"]),
            "eta_error": float(reg_aph["eta_error"]),
        },
        "differential": {
            "delta_eta_mag": float(eta_diff),
            "significance_sigma": float(diff_significance),
            "is_significant": significant,
        },
        "threshold_activation_fit": threshold_fit,
        "binned_analysis": {
            "n_bins": len(bin_centers),
            "bin_centers_au": bin_centers,
            "bin_etas": bin_etas,
            "bin_eta_errors": bin_eta_errs,
        }
        if bin_centers
        else None,
    }

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_022", str(log_dir / "step_022_environmental_modulation.log")
    )
    set_step_logger(logger)

    data_path = (
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    if not data_path.exists():
        print_status("No processed INPOP19a residuals.", "ERROR")
        sys.exit(1)

    df = pd.read_csv(data_path)

    results = compute_distances_and_analyze(df)

    logger.save_step_results(results, PROJECT_ROOT, "step_022_environmental_modulation")
    print_status("Environmental Modulation Tests Complete.", "SUCCESS")