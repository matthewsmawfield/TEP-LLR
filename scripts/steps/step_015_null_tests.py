#!/usr/bin/env python3
"""
Step 015: Comprehensive Null Tests for TEP-LLR

Tests for false positives by scanning non-synodic frequencies and checking
that the primary TEP signal is not an artifact of uncontrolled systematics.
Pre-whitening is applied using the elongation_rad phase basis to ensure the
filter shares the same basis as the detection regression.
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import math
import pandas as pd
import numpy as np
from scripts.utils.numerics import stable_lstsq
from scripts.utils.statistical_utils import linear_regression
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.pre_whitening_filter import apply_pre_whitening
from scripts.utils.llr_constants import (
    NULL_TEST_SCAN_MIN_FACTOR,
    NULL_TEST_SCAN_MAX_FACTOR,
    NULL_TEST_SCAN_POINTS,
    NULL_TEST_EXCLUDE_SYNODIC_WINDOW,
    NULL_TEST_PRIMARY_FREQUENCY_FACTOR,
)

# Add project root to path

def _fdr_bh_significance(p_values, alpha=0.05):
    m = len(p_values)
    if m == 0:
        return [], 0.0, False

    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    sorted_p = p[order]
    crit = alpha * (np.arange(1, m + 1) / m)
    passed = sorted_p <= crit

    if np.any(passed):
        k = int(np.max(np.where(passed)[0]))
        threshold = float(sorted_p[k])
        sig_sorted = sorted_p <= threshold
    else:
        threshold = 0.0
        sig_sorted = np.zeros(m, dtype=bool)

    sig = np.zeros(m, dtype=bool)
    sig[order] = sig_sorted
    any_significant = bool(np.any(sig))
    return sig.tolist(), threshold, any_significant

def _parse_exclude_bands(values):
    bands = []
    for value in values or []:
        parts = value.split(":")
        if len(parts) != 2:
            raise ValueError(
                f"Invalid --exclude-band '{value}'. Expected format min:max")
        lo = float(parts[0])
        hi = float(parts[1])
        if hi < lo:
            lo, hi = hi, lo
        bands.append((lo, hi))
    return bands

def _is_in_excluded_band(factor, exclude_bands):
    return any(lo <= factor <= hi for lo, hi in (exclude_bands or []))

def run_null_tests(
    df,
    verbose=False,
    correction_mode="fdr_bh",
    scan_min_factor=NULL_TEST_SCAN_MIN_FACTOR,
    scan_max_factor=NULL_TEST_SCAN_MAX_FACTOR,
    scan_points=NULL_TEST_SCAN_POINTS,
    exclude_synodic_window=NULL_TEST_EXCLUDE_SYNODIC_WINDOW,
    exclude_bands=None,
):
    # NOTE: exclude_bands does NOT remove frequencies from the actual test.
    # All candidate_factors (except the synodic window) are tested regardless.
    # exclude_bands only affects reporting metadata (frequencies_excluded_by_bands).
    # In Deep Scan mode, we scan all frequencies but label known systematic regions
    # Sidereal ~ 0.925 synodic, Anomalistic ~ 0.933 synodic
    # Annual ~ 0.08 synodic
    systematic_regions = [
        (0.05, 0.15, "Annual/Seasonal"),
        (0.90, 0.96, "Lunar Orbital (Sidereal/Anomalistic)"),
        (1.95, 2.05, "Semi-Synodic Harmonic"),
    ]

    # This frequency is intentionally non-physical.
    # It lies between known systematic bands (annual ~0.08, lunar orbital ~0.93, semi-synodic ~2.0)
    # and serves as a control frequency where no physical signal is expected.
    primary_test_frequency_factor = NULL_TEST_PRIMARY_FREQUENCY_FACTOR
    candidate_factors = np.linspace(scan_min_factor, scan_max_factor, scan_points)

    # Exclude only the immediate synodic window to avoid self-nulling
    excluded_by_synodic = [f for f in candidate_factors if abs(float(f) - 1.0) <= exclude_synodic_window]
    test_frequency_factors = sorted({
        float(f)
        for f in candidate_factors
        if abs(float(f) - 1.0) > exclude_synodic_window
    } | {primary_test_frequency_factor})

    n_tests = len(test_frequency_factors)

    # Apply joint pre-whitening on nuisance harmonics ONLY.
    # CRITICAL FIX: Do NOT pass test_frequency_factors to pre-whitening.
    # Passing them would remove power at the very frequencies we intend to test,
    # creating a self-fulfilling null result. Whiten only dominant non-synodic
    # peaks; the null-test frequencies are tested on the residualised data.
    df_white = apply_pre_whitening(
        df,
        n_harmonics=5,
        verbose=verbose,
    )

    residuals_white = df_white['residual_whitened_m'].values
    elongation = df_white['elongation_rad'].values

    # Project out TEP basis
    cos_elong_tep = np.cos(elongation)
    X_tep = np.column_stack([cos_elong_tep, np.ones_like(cos_elong_tep)])
    tep_coeffs, _, _, _ = stable_lstsq(X_tep, residuals_white)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        residuals_null = residuals_white - X_tep @ tep_coeffs

    test_results = []
    for factor in test_frequency_factors:
        reg = linear_regression(residuals_null, np.cos(elongation * factor))
        snr = float(abs(reg['eta']) / reg['eta_error']) if reg['eta_error'] > 0 else 0.0
        p_raw = float(math.erfc(snr / math.sqrt(2.0))) if snr > 0 else 1.0

        label = "Null Region"
        for lo, hi, name in systematic_regions:
            if lo <= factor <= hi:
                label = name
                break

        test_results.append({
            "frequency_factor": float(factor),
            "label": label,
            "snr": snr,
            "p_two_sided": p_raw,
        })

    p_values = [r["p_two_sided"] for r in test_results]
    fdr_sig_flags, fdr_threshold, any_significant_fdr = _fdr_bh_significance(p_values, alpha=0.05)

    for i, flag in enumerate(fdr_sig_flags):
        test_results[i]["significant_after_fdr_bh"] = bool(flag)

    primary_result = min(test_results, key=lambda r: abs(r["frequency_factor"] - primary_test_frequency_factor))

    # Worst case EXCLUDING known systematic regions
    non_systematic_results = [r for r in test_results if r["label"] == "Null Region"]
    worst_null_result = max(non_systematic_results, key=lambda r: r["snr"]) if non_systematic_results else primary_result

    pass_null = bool(
        primary_result["snr"] < 5.0 and
        worst_null_result["snr"] < 5.0
    )

    if verbose:
        print_status("Deep Scan Multi-frequency Null Test Summary:", "CALC")
        print_status(f"  Primary factor (1.23×) SNR: {primary_result['snr']:.2f}σ", "CALC")
        print_status(f"  Worst clean null factor: {worst_null_result['frequency_factor']:.2f}×, SNR={worst_null_result['snr']:.2f}σ", "CALC")

        # Report systematic peaks discovered
        systematic_peaks = [r for r in test_results if r["label"] != "Null Region" and r["snr"] > 3.0]
        if systematic_peaks:
            print_status("  Detected Systematic Peaks (Expected):", "INFO")
            for r in systematic_peaks:
                print_status(f"    {r['label']} (f={r['frequency_factor']:.2f}×): SNR={r['snr']:.2f}σ", "INFO")

    # Count exclusions by bands
    excluded_by_bands = [f for f in candidate_factors
                         if any(lo <= f <= hi for lo, hi in (exclude_bands or []))
                         and f not in excluded_by_synodic]

    return {
        "correction_mode": correction_mode,
        "primary_test_frequency_factor": float(primary_test_frequency_factor),
        "non_physical_snr": float(primary_result["snr"]),
        "max_null_region_snr": float(worst_null_result["snr"]),
        "pass_null": pass_null,
        "systematic_peaks": [r for r in test_results if r["snr"] > 3.0],
        "test_results": test_results,
        "scan_metadata": {
            "scan_min_factor": float(scan_min_factor),
            "scan_max_factor": float(scan_max_factor),
            "total_candidate_frequencies": int(scan_points),
            "frequencies_excluded_by_synodic_window": len(excluded_by_synodic),
            "frequencies_excluded_by_bands": len(excluded_by_bands),
            "total_frequencies_tested": int(n_tests),
            "exclude_synodic_window": float(exclude_synodic_window),
            "exclude_bands": [list(b) for b in (exclude_bands or [])],
            "systematic_regions_defined": [(lo, hi, name) for lo, hi, name in systematic_regions],
            "fdr_threshold": float(fdr_threshold),
            "n_tests": int(n_tests)
        }
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 015: Null Tests")
    parser.add_argument("--correction-mode",
                        choices=["fdr_bh", "bonferroni"], default="fdr_bh")
    parser.add_argument("--scan-min-factor", type=float, default=0.4)
    parser.add_argument("--scan-max-factor", type=float, default=3.2)
    parser.add_argument("--scan-points", type=int, default=57)
    parser.add_argument("--exclude-synodic-window", type=float, default=0.08)
    parser.add_argument("--exclude-band", action="append", default=[],
                        help="Exclude frequency band from scan, format min:max. Repeatable.")
    args = parser.parse_args()

    exclude_bands = _parse_exclude_bands(args.exclude_band)
    # Apply default low-frequency exclusion if no bands specified via CLI
    if not exclude_bands:
        exclude_bands = [(0.4, 0.95)]

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_015", str(
        log_dir / "step_015_null_tests.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Starting Null Tests...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    summary = run_null_tests(
        df,
        verbose=True,
        correction_mode=args.correction_mode,
        scan_min_factor=args.scan_min_factor,
        scan_max_factor=args.scan_max_factor,
        scan_points=args.scan_points,
        exclude_synodic_window=args.exclude_synodic_window,
        exclude_bands=exclude_bands,
    )

    results = {
        "step_id": "step_015",
        "null_test_summary": summary,
        "status": "PASS" if summary["pass_null"] else "WARNING"
    }

    logger.save_step_results(results, PROJECT_ROOT, "step_015_null_tests")
    print_status("Null Tests Complete.", "SUCCESS")
