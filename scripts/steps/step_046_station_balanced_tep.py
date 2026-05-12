#!/usr/bin/env python3
"""
Step 046: Station-Balanced TEP Analysis

Addresses the Grasse-dominance concern (74% of observations) by creating
subsamples where stations contribute equally, then testing whether the TEP
signal persists.  This is the most direct falsification test: if the signal
is driven by Grasse-specific systematics, it should vanish when Grasse is
downweighted to match other stations.

Methods:
  1. Equal-N subsample: Each station contributes min(N_stations) observations
  2. Grasse-capped: Grasse downsampled to match the next-largest station (APO)
  3. Stratified bootstrap: Multiple iterations of balanced subsampling
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import argparse
import numpy as np
import pandas as pd
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status


def create_equal_n_subsample(df, min_obs_per_station=None, seed=42):
    """Create a balanced subsample with equal observations per station."""
    rng = np.random.RandomState(seed)
    stations = df['station'].unique()

    if min_obs_per_station is None:
        min_obs_per_station = min(len(df[df['station'] == s]) for s in stations)

    balanced = []
    for s in stations:
        sdf = df[df['station'] == s]
        if len(sdf) >= min_obs_per_station:
            idx = rng.choice(len(sdf), size=min_obs_per_station, replace=False)
            balanced.append(sdf.iloc[idx])
        else:
            balanced.append(sdf)

    return pd.concat(balanced, ignore_index=True)


def create_grasse_capped_subsample(df, cap_station='APO', seed=42):
    """Downsample Grasse to match the cap_station's observation count."""
    rng = np.random.RandomState(seed)
    cap_n = len(df[df['station'] == cap_station])

    balanced = []
    for s in df['station'].unique():
        sdf = df[df['station'] == s]
        if s == 'Grasse' and len(sdf) > cap_n:
            idx = rng.choice(len(sdf), size=cap_n, replace=False)
            balanced.append(sdf.iloc[idx])
        else:
            balanced.append(sdf)

    return pd.concat(balanced, ignore_index=True)


def run_balanced_analysis(df, name, outlier_threshold=6.0):
    """Run TEP regression on a balanced subsample."""
    residuals = df['residual_m'].values
    cos_elong = np.cos(df['elongation_rad'].values)

    # Outlier cleaning
    outlier_mask = detect_outliers_sigma(residuals, sigma_threshold=outlier_threshold)
    kept = ~outlier_mask
    n_outliers = int(outlier_mask.sum())

    if kept.sum() < 100:
        return None

    reg = linear_regression(residuals[kept], cos_elong[kept])
    snr = abs(reg['eta']) / reg['eta_error'] if reg['eta_error'] > 0 else 0.0

    return {
        'name': name,
        'n_total': len(df),
        'n_used': int(kept.sum()),
        'n_outliers': n_outliers,
        'eta': float(reg['eta']),
        'eta_error': float(reg['eta_error']),
        'snr': float(snr),
        'significant': bool(snr >= 3.0),
        'stations': {s: int((df['station'] == s).sum()) for s in df['station'].unique()}
    }


def main():
    parser = argparse.ArgumentParser(description="Step 046: Station-Balanced TEP Analysis")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_046", str(log_dir / "step_046_station_balanced_tep.log"))
    set_step_logger(logger)

    print_status("Starting Step 046: Station-Balanced TEP Analysis...", "TITLE")
    print_status("Purpose: Test TEP signal persistence when Grasse dominance is removed", "INFO")

    input_path = PROJECT_ROOT / 'data' / 'processed' / 'INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    print_status(f"Loaded {len(df):,} observations", "INFO")
    for s in sorted(df['station'].unique()):
        print_status(f"  Station {s}: {(df['station'] == s).sum():,} obs", "INFO")

    results = []

    # 1. Full sample baseline
    print_status("", "INFO")
    print_status(">>> TEST 1: Full-sample baseline", "PROCESS")
    full = run_balanced_analysis(df, 'full_sample')
    if full:
        print_status(f"    η = {full['eta']:.3e} ± {full['eta_error']:.3e} ({full['snr']:.2f}σ)", "CALC")
        results.append(full)

    # 2. Equal-N subsample
    print_status("", "INFO")
    print_status(">>> TEST 2: Equal-N subsample (all stations = min N)", "PROCESS")
    eq_df = create_equal_n_subsample(df, seed=42)
    eq = run_balanced_analysis(eq_df, 'equal_n_subsample')
    if eq:
        print_status(f"    η = {eq['eta']:.3e} ± {eq['eta_error']:.3e} ({eq['snr']:.2f}σ)", "CALC")
        for s, n in sorted(eq['stations'].items()):
            print_status(f"      {s}: {n} obs", "INFO")
        results.append(eq)

    # 3. Grasse-capped subsample
    print_status("", "INFO")
    print_status(">>> TEST 3: Grasse-capped subsample (Grasse = APO N)", "PROCESS")
    cap_df = create_grasse_capped_subsample(df, cap_station='APO', seed=42)
    cap = run_balanced_analysis(cap_df, 'grasse_capped')
    if cap:
        print_status(f"    η = {cap['eta']:.3e} ± {cap['eta_error']:.3e} ({cap['snr']:.2f}σ)", "CALC")
        for s, n in sorted(cap['stations'].items()):
            print_status(f"      {s}: {n} obs", "INFO")
        results.append(cap)

    # 4. Stratified bootstrap: multiple equal-N iterations
    print_status("", "INFO")
    print_status(">>> TEST 4: Stratified bootstrap (200 iterations of equal-N)", "PROCESS")
    bootstrap_snrs = []
    bootstrap_etas = []
    N_BOOT = 200
    for i in range(N_BOOT):
        b_df = create_equal_n_subsample(df, seed=42 + i)
        b = run_balanced_analysis(b_df, f'bootstrap_{i}')
        if b:
            bootstrap_snrs.append(b['snr'])
            bootstrap_etas.append(b['eta'])

    if bootstrap_etas:
        eta_mean = np.mean(bootstrap_etas)
        eta_std = np.std(bootstrap_etas)
        snr_mean = np.mean(bootstrap_snrs)
        # 95% CI using percentile method
        ci_lower = np.percentile(bootstrap_etas, 2.5)
        ci_upper = np.percentile(bootstrap_etas, 97.5)
        print_status(f"    Mean η = {eta_mean:.3e} ± {eta_std:.3e}", "CALC")
        print_status(f"    95% CI = [{ci_lower:.3e}, {ci_upper:.3e}]", "CALC")
        print_status(f"    Mean SNR = {snr_mean:.2f}σ", "CALC")
        print_status(f"    All significant? {all(s >= 3.0 for s in bootstrap_snrs)}", "CALC")

    # Summary
    print_status("", "INFO")
    print_status("═══ SUMMARY", "TITLE")
    all_significant = all(r['significant'] for r in results)
    if all_significant:
        print_status("All balanced subsamples detect TEP at >= 3σ", "SUCCESS")
        print_status("Conclusion: Signal is NOT driven by Grasse dominance", "SUCCESS")
    else:
        print_status("Some balanced subsamples fail to detect TEP", "WARNING")
        print_status("Conclusion: Signal may be sensitive to station weighting", "WARNING")

    output = {
        "step_id": "step_046",
        "status": "PASS" if all_significant else "WARNING",
        "tests": results,
        "bootstrap": {
            "n_iterations": N_BOOT,
            "eta_mean": float(np.mean(bootstrap_etas)) if bootstrap_etas else None,
            "eta_std": float(np.std(bootstrap_etas)) if bootstrap_etas else None,
            "eta_ci95_lower": float(np.percentile(bootstrap_etas, 2.5)) if bootstrap_etas else None,
            "eta_ci95_upper": float(np.percentile(bootstrap_etas, 97.5)) if bootstrap_etas else None,
            "snr_mean": float(np.mean(bootstrap_snrs)) if bootstrap_snrs else None,
            "snr_std": float(np.std(bootstrap_snrs)) if bootstrap_snrs else None,
            "all_significant": all(s >= 3.0 for s in bootstrap_snrs) if bootstrap_snrs else None
        },
        "conclusion": "Signal persists in station-balanced subsamples" if all_significant else "Signal weakens when station balance is enforced"
    }

    logger.save_step_results(output, PROJECT_ROOT, "step_046_station_balanced_tep")
    print_status("Station-Balanced Analysis Complete.", "SUCCESS")


if __name__ == "__main__":
    main()
