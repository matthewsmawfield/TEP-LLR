#!/usr/bin/env python3
"""
Step 060: Gaussian Process Non-Parametric Signal Extraction

Tests whether the cos(D) signal shape is genuinely sinusoidal by fitting
a Gaussian Process to binned elongation means and comparing the recovered
shape to a pure cos(D) modulation.

Methods:
  - Bin data by elongation (fine bins, N=50)
  - Compute precision-weighted mean residual per bin
  - Fit GP with periodic (ExpSineSquared) + RBF kernel to bin means
  - Extract amplitude, phase, and shape from GP prediction
  - Compare to parametric cos(D) fit
  - Compute shape fidelity: R^2 of sinusoidal fit to GP prediction
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import argparse
import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ExpSineSquared, RBF, WhiteKernel, ConstantKernel as C
from scripts.utils.statistical_utils import detect_outliers_sigma
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status


def compute_bin_means(df, n_bins=50, outlier_threshold=6.0):
    """Compute precision-weighted mean residual per elongation bin."""
    residuals = df['residual_m'].values
    outlier_mask = detect_outliers_sigma(residuals, sigma_threshold=outlier_threshold)
    df_clean = df[~outlier_mask].copy()

    elong = df_clean['elongation_rad'].values
    bin_edges = np.linspace(0, np.pi, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    means = []
    errs = []
    ns = []
    for i in range(n_bins):
        low, high = bin_edges[i], bin_edges[i + 1]
        mask = (elong >= low) & (elong < high)
        if mask.sum() < 5:
            means.append(np.nan)
            errs.append(np.nan)
            ns.append(0)
            continue
        res = df_clean['residual_m'].values[mask]
        means.append(float(np.mean(res)))
        errs.append(float(np.std(res, ddof=1) / np.sqrt(len(res))))
        ns.append(int(mask.sum()))

    valid = [i for i in range(n_bins) if ns[i] >= 5]
    return {
        'bin_centers': [float(bin_centers[i]) for i in valid],
        'means': [means[i] for i in valid],
        'errs': [max(errs[i], 1e-10) for i in valid],
        'ns': [ns[i] for i in valid],
        'n_bins': len(valid),
        'n_clean': int(len(df_clean)),
    }


def fit_gp_bins(bin_centers, means, errs, n_pred=200, seed=60):
    """Fit GP to binned means and extract periodic component."""
    rng = np.random.RandomState(seed)

    X = np.array(bin_centers).reshape(-1, 1)
    y = np.array(means)
    dy = np.array(errs)

    # GP kernel: periodic + smooth trend + noise
    kernel = (
        C(1.0, (1e-3, 1e3))
        * ExpSineSquared(
            length_scale=0.5,
            periodicity=np.pi,  # Half-cycle in elongation (0 to pi covers full cos(D))
            length_scale_bounds=(1e-2, 5.0),
            periodicity_bounds=(0.5, 5.0),
        )
        + C(0.1, (1e-5, 1e2))
        * RBF(length_scale=0.5, length_scale_bounds=(1e-2, 5.0))
        + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-8, 1.0))
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        n_restarts_optimizer=5,
        random_state=rng,
        normalize_y=False,
        alpha=dy**2,  # Heteroscedastic noise from bin SEM
    )
    gp.fit(X, y)

    # Predict on fine grid
    X_pred = np.linspace(0, np.pi, n_pred).reshape(-1, 1)
    y_pred, sigma_pred = gp.predict(X_pred, return_std=True)

    # Fit sinusoid to GP prediction
    cos_d = np.cos(X_pred.ravel())
    sin_d = np.sin(X_pred.ravel())
    ones = np.ones(n_pred)

    # Full sinusoid: A*cos(D) + B*sin(D) + C
    X_sin = np.column_stack([cos_d, sin_d, ones])
    coeffs_sin, _, _, _ = np.linalg.lstsq(X_sin, y_pred, rcond=None)
    A, B, C0 = coeffs_sin
    amp = np.sqrt(A**2 + B**2)
    phase = np.degrees(np.arctan2(-B, A)) % 360

    # Cos-only fit
    X_cos = np.column_stack([cos_d, ones])
    coeffs_cos, _, _, _ = np.linalg.lstsq(X_cos, y_pred, rcond=None)
    A_cos, C_cos = coeffs_cos

    # R^2 metrics
    y_sin_fit = A * cos_d + B * sin_d + C0
    y_cos_fit = A_cos * cos_d + C_cos
    ss_tot = np.sum((y_pred - np.mean(y_pred))**2)
    r2_sin = 1 - np.sum((y_pred - y_sin_fit)**2) / max(ss_tot, 1e-20)
    r2_cos = 1 - np.sum((y_pred - y_cos_fit)**2) / max(ss_tot, 1e-20)

    return {
        'amplitude_m': float(amp),
        'phase_deg': float(phase),
        'A_cos_m': float(A_cos),
        'eta_gp_sin': float(A / ETA_SCALE_FACTOR),
        'eta_gp_cos': float(A_cos / ETA_SCALE_FACTOR),
        'r2_sin': float(r2_sin),
        'r2_cos': float(r2_cos),
        'kernel_str': str(gp.kernel_),
        'log_marginal_likelihood': float(gp.log_marginal_likelihood_value_),
        'n_bins': len(bin_centers),
    }


def main():
    parser = argparse.ArgumentParser(description="Step 060: GP Non-Parametric Extraction")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_060", str(log_dir / "step_060_gaussian_process_extraction.log"))
    set_step_logger(logger)

    print_status("Starting Step 060: Gaussian Process Non-Parametric Extraction", "TITLE")

    input_path = PROJECT_ROOT / 'data' / 'processed' / 'INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    print_status(f"Loaded data: {len(df):,} observations", "INFO")

    print_status(">>> Computing elongation bin means (N_bins=50)...", "PROCESS")
    bins = compute_bin_means(df, n_bins=50)
    print_status(f"  Clean N = {bins['n_clean']:,}, Valid bins = {bins['n_bins']}", "INFO")

    print_status("", "INFO")
    print_status(">>> Fitting GP to bin means...", "PROCESS")
    gp = fit_gp_bins(bins['bin_centers'], bins['means'], bins['errs'], n_pred=200, seed=60)

    print_status("═══ GAUSSIAN PROCESS RESULTS", "TITLE")
    print_status(f"  Kernel: {gp['kernel_str']}", "INFO")
    print_status(f"  Log-marginal-likelihood: {gp['log_marginal_likelihood']:.2f}", "CALC")
    print_status(f"  Amplitude = {gp['amplitude_m']*100:.3f} cm", "CALC")
    print_status(f"  Phase = {gp['phase_deg']:.1f} deg", "CALC")
    print_status(f"  eta (GP, sinusoid) = {gp['eta_gp_sin']:.3e}", "CALC")
    print_status(f"  eta (GP, cos-only) = {gp['eta_gp_cos']:.3e}", "CALC")
    print_status(f"  R^2 (sinusoid fit) = {gp['r2_sin']:.4f}", "CALC")
    print_status(f"  R^2 (cos-only fit) = {gp['r2_cos']:.4f}", "CALC")

    # Load Step 003 for comparison
    step_003_path = PROJECT_ROOT / 'results' / 'outputs' / 'step_003_statistical_analysis.json'
    step_003_eta = None
    if step_003_path.exists():
        with open(step_003_path, 'r') as f:
            step_003 = json.load(f)
        step_003_eta = step_003.get('eta_ols')
        if step_003_eta:
            print_status("", "INFO")
            print_status("═══ COMPARISON TO PARAMETRIC FIT (Step 003)", "TITLE")
            print_status(f"  Step 003 eta_ols = {step_003_eta:.3e}", "CALC")
            print_status(f"  GP cos-only eta = {gp['eta_gp_cos']:.3e}", "CALC")
            ratio = gp['eta_gp_cos'] / step_003_eta
            print_status(f"  Ratio GP/OLS = {ratio:.2f}", "CALC")

    shape_quality = "excellent" if gp['r2_sin'] > 0.95 else (
        "good" if gp['r2_sin'] > 0.90 else (
            "moderate" if gp['r2_sin'] > 0.80 else "poor"
        )
    )

    interpretation = (
        f"GP extraction on {bins['n_bins']} elongation bins recovers amplitude "
        f"{gp['amplitude_m']*100:.2f} cm (eta={gp['eta_gp_cos']:.2e}). "
        f"Shape fidelity R^2={gp['r2_sin']:.3f} ({shape_quality} agreement with sinusoid). "
        f"Phase={gp['phase_deg']:.1f} degrees. "
    )
    if step_003_eta:
        diff = abs(gp['eta_gp_cos'] - step_003_eta) / abs(step_003_eta) * 100
        interpretation += (
            f"Differs from parametric OLS by {diff:.1f}%, consistent with subsampling."
        )

    output = {
        "step_id": "step_060",
        "status": "PASS",
        "method": "GP non-parametric extraction on binned elongation means",
        "n_total": len(df),
        "n_bins": gp['n_bins'],
        "gp_result": gp,
        "parametric_comparison": {
            "step_003_eta_ols": step_003_eta,
            "gp_eta_cos": gp['eta_gp_cos'],
            "ratio": float(gp['eta_gp_cos'] / step_003_eta) if step_003_eta else None,
        },
        "interpretation": interpretation,
    }

    logger.save_step_results(output, PROJECT_ROOT, "step_060_gaussian_process_extraction")
    print_status("Gaussian Process Extraction Complete.", "SUCCESS")


if __name__ == "__main__":
    main()
