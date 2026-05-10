#!/usr/bin/env python3
"""
Step 012: Noise/Signal Injection Tests for TEP-LLR
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import argparse
import pandas as pd
import numpy as np
from scripts.utils.statistical_utils import linear_regression
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

# Add the project root to the Python path

def run_injection_test(df, verbose=False):
    print_status("="*60, "INFO")
    print_status("NOISE/SIGNAL INJECTION TESTS - PIPELINE VALIDATION", "TITLE")
    print_status("WARNING: THESE RESULTS ARE BASED ON SYNTHETIC DATA", "WARNING")
    print_status("PURPOSE: ALGORITHM ROBUSTNESS TESTING ONLY", "INFO")
    print_status("="*60, "INFO")

    residuals = df['residual_m'].values
    elongation = df['elongation_rad'].values
    cos_elong = np.cos(elongation)
    n = len(residuals)

    rms_original = np.std(residuals)

    print_status(f"[DATA] Original dataset: N={n:,} observations", "INFO")
    print_status(f"[DATA] Original RMS: {rms_original:.6f} m", "INFO")
    print_status(f"[DATA] Original residual range: [{np.min(residuals):.6f}, {np.max(residuals):.6f}] m", "INFO")

    # TEST 1: Null test - Shuffle residuals (destroys any real correlation)
    print_status("", "INFO")
    print_status("TEST 1: NULL TEST (SHUFFLED RESIDUALS)", "PROCESS")
    print_status("  Purpose: Verify that destroying phase correlation eliminates signal", "INFO")

    np.random.seed(42)  # Reproducibility
    res_shuffled = np.random.permutation(residuals)
    reg_null = linear_regression(res_shuffled, cos_elong)
    snr_null = abs(reg_null['eta']) / reg_null['eta_error']

    print_status("  [CALC] Linear regression on shuffled data:", "CALC")
    print_status(f"  [CALC]    η estimate: {reg_null['eta']:.6e}", "CALC")
    print_status(f"  [CALC]    SNR: {snr_null:.2f}σ", "CALC")
    print_status("  [CALC]    Expected: SNR < 3 (no signal in shuffled data)", "CALC")
    print_status(f"  [RESULT] {'PASS ✓' if snr_null < 3 else 'FAIL ✗'}", "SUCCESS" if snr_null < 3 else "WARNING")

    # TEST 2: Add Gaussian noise at various levels
    print_status("", "INFO")
    print_status("TEST 2: NOISE ROBUSTNESS (ADDITIVE GAUSSIAN NOISE)", "PROCESS")
    print_status("  Purpose: Determine how much noise signal can survive", "INFO")

    noise_levels = [0.5, 1.0, 2.0, 3.0]  # multiples of RMS (standard test range: 0.5-3x baseline noise)
    noise_results = []

    for noise_mult in noise_levels:
        np.random.seed(42 + int(noise_mult * 10))  # Unique seed per level
        noise = np.random.normal(0, noise_mult * rms_original, n)
        res_noisy = residuals + noise
        reg_noisy = linear_regression(res_noisy, cos_elong)
        snr_noisy = abs(reg_noisy['eta']) / reg_noisy['eta_error']

        noise_results.append({
            "noise_multiplier": noise_mult,
            "total_rms": float(np.std(res_noisy)),
            "eta": float(reg_noisy['eta']),
            "eta_error": float(reg_noisy['eta_error']),
            "snr": float(snr_noisy),
            "significant": bool(snr_noisy > 3)
        })

        print_status(f"  [CALC] Noise level {noise_mult}× RMS: η = {reg_noisy['eta']:.6e}, SNR = {snr_noisy:.2f}σ {'✓' if snr_noisy > 3 else ''}", "SUCCESS" if snr_noisy > 3 else "CALC")

    # TEST 3: Signal recovery test (inject known signal into pure noise)
    print_status("", "INFO")
    print_status("TEST 3: SIGNAL RECOVERY (INJECTION INTO PURE NOISE)", "PROCESS")
    print_status("  Purpose: Validate pipeline can recover known injected signals", "INFO")

    np.random.seed(42)
    pure_noise = np.random.normal(0, rms_original, n)

    # Load measured eta from step_002 output for injection test
    step_002_path = PROJECT_ROOT / 'results' / 'outputs' / 'step_002_statistical_analysis.json'
    if step_002_path.exists():
        with open(step_002_path, 'r') as f:
            step_002_results = json.load(f)
        eta_injected = step_002_results.get('eta_ols', 0)
        print_status(f"  [DATA] Loaded measured η from step_002: {eta_injected:.4e}", "INFO")
    else:
        raise FileNotFoundError(f"Step 002 results not found: {step_002_path}. Run pipeline step 002 first.")

    signal_injected = 13.0 * eta_injected * cos_elong
    res_with_signal = pure_noise + signal_injected

    reg_recovery = linear_regression(res_with_signal, cos_elong)
    eta_recovered = reg_recovery['eta']
    recovery_error = abs(eta_recovered - eta_injected) / \
        abs(eta_injected) * 100

    print_status("  [CALC] Injected signal parameters:", "CALC")
    print_status(f"  [CALC]    η_injected: {eta_injected:.6e}", "CALC")
    print_status("  [CALC] Recovery results:", "CALC")
    print_status(f"  [CALC]    η_recovered: {eta_recovered:.6e}", "CALC")
    print_status(f"  [CALC]    Recovery error: {recovery_error:.1f}%", "CALC")
    print_status(f"  [CALC]    SNR: {abs(eta_recovered)/reg_recovery['eta_error']:.2f}σ", "CALC")
    print_status(f"  [RESULT] {'PASS ✓' if recovery_error < 20 else 'WARNING'} (error < 20%)", "SUCCESS" if recovery_error < 20 else "WARNING")

    # TEST 4: Detection threshold analysis
    print_status("", "INFO")
    print_status("TEST 4: DETECTION THRESHOLD ANALYSIS", "PROCESS")
    print_status("  Purpose: Compute minimum detectable η at various confidence levels", "INFO")

    confidence_levels = [1.0, 2.0, 3.0, 5.0]  # sigma
    threshold_results = []

    for sigma in confidence_levels:
        # Minimum detectable eta: SNR = sigma
        # sigma = |eta| / eta_error => |eta| = sigma * eta_error
        # eta_error depends on sample size and noise
        # eta_error ≈ sigma_resid / (sqrt(N) * 13 * rms_cos_elong)
        rms_cos = np.std(cos_elong)
        eta_min = sigma * rms_original / (np.sqrt(n) * 13.0 * rms_cos)

        threshold_results.append({
            "confidence_sigma": sigma,
            "min_detectable_eta": float(eta_min)
        })

        if verbose:
            print_status(f"  [CALC] {sigma}σ detection threshold:", "CALC")
            print_status(f"  [CALC]    min |η| = {eta_min:.6e}", "CALC")
            detected = abs(reg_recovery['eta']) > eta_min
            print_status(
                f"  [CALC]    Current signal detected: {'YES ✓' if detected else 'NO'}", "SUCCESS" if detected else "WARNING")

    if verbose:
        print_status("", "INFO")
        print_status("="*60, "INFO")
        print_status("NOISE INJECTION SUMMARY", "TITLE")
        print_status("="*60, "INFO")
        print_status(
            f"  Null test:           SNR = {snr_null:.2f}σ (should be < 3) {'✓' if snr_null < 3 else '✗'}", "SUCCESS" if snr_null < 3 else "WARNING")
        print_status(f"  Noise at 2× RMS:     Signal {'survives' if any(r['noise_multiplier'] == 2.0 and r['significant'] for r in noise_results) else 'lost'} {'✓' if any(r['noise_multiplier'] == 2.0 and r['significant'] for r in noise_results) else ''}", "SUCCESS" if any(
            r['noise_multiplier'] == 2.0 and r['significant'] for r in noise_results) else "WARNING")
        print_status(f"  Signal recovery:     {recovery_error:.1f}% error {'✓' if recovery_error < 20 else '✗'}",
                     "SUCCESS" if recovery_error < 20 else "WARNING")
        print_status("="*60, "INFO")

    return {
        "null_test_eta": float(reg_null['eta']),
        "null_test_eta_error": float(reg_null['eta_error']),
        "null_test_snr": float(snr_null),
        "null_test_pass": bool(snr_null < 3),
        "noise_robustness_tests": noise_results,
        "signal_recovery": {
            "eta_injected": float(eta_injected),
            "eta_recovered": float(eta_recovered),
            "recovery_error_percent": float(recovery_error),
            "recovery_pass": bool(recovery_error < 20)
        },
        "detection_thresholds": threshold_results,
        "injection_valid": bool(snr_null < 3 and recovery_error < 20)
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 011: Noise Signal Injection")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_011", str(
        log_dir / "step_011_noise_signal_injection.log"))
    set_step_logger(logger)

    print_status("Starting Noise/Signal Injection Tests...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    summary = run_injection_test(df)

    results = {
        "step_id": "step_011",
        "data_type": "SYNTHETIC (PIPELINE VALIDATION)",
        "injection_results": summary,
        "status": "PASS" if summary["injection_valid"] else "WARNING"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_011_noise_signal_injection")
    print_status("Injection Tests Complete.", "SUCCESS")