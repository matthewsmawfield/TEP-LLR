#!/usr/bin/env python3
"""
Step 052: Station Sample and Distribution Analysis
===================================================
Properly analyze station sample sizes, elongation coverage, temporal coverage,
and joint distributions to understand prediction failure and signal stability.

This addresses the critique that negative predictive R^2 in cross-validation
may arise from covariate shift (different stations/elongations in train vs test)
rather than signal non-existence.
"""

import sys
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
import numpy as np
from scripts.utils.numerics import stable_lstsq
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning
import pandas as pd
from scipy import stats

log_dir = PROJECT_ROOT / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_052", str(log_dir / "step_052_station_distribution_analysis.log"))
set_step_logger(logger)


def load_data():
    df = pd.read_csv(PROJECT_ROOT / "data/processed/INPOP19a_all_stations_residuals.csv")
    # Recompute cosD, annual, monthly terms
    df['cosD'] = np.cos(df['elongation_rad'].values)
    df['cos2d'] = np.cos(2 * df['elongation_rad'].values)
    year = df['date_julian'].values / 365.25
    df['sin_y'] = np.sin(2 * np.pi * year)
    df['cos_y'] = np.cos(2 * np.pi * year)
    month = df['date_julian'].values / 27.32
    df['sin_m'] = np.sin(2 * np.pi * month)
    df['cos_m'] = np.cos(2 * np.pi * month)
    return df


def analyze_station_coverage(df):
    """Analyze per-station sample sizes, date ranges, and phase coverage."""
    stations = df['station'].unique()
    results = {}
    for stn in stations:
        sub = df[df['station'] == stn]
        n = len(sub)
        date_min = float(sub['date_julian_year'].min())
        date_max = float(sub['date_julian_year'].max())
        elong_mean = float(sub['elongation_rad'].mean())
        elong_std = float(sub['elongation_rad'].std())
        cosD_mean = float(sub['cosD'].mean())
        cosD_std = float(sub['cosD'].std())
        rms = float(sub['residual_m'].std())
        # Phase coverage: fraction in each quadrant
        q1 = float(np.mean((sub['elongation_rad'] >= 0) & (sub['elongation_rad'] < np.pi/2)))
        q2 = float(np.mean((sub['elongation_rad'] >= np.pi/2) & (sub['elongation_rad'] < np.pi)))
        q3 = float(np.mean((sub['elongation_rad'] >= np.pi) & (sub['elongation_rad'] < 3*np.pi/2)))
        q4 = float(np.mean((sub['elongation_rad'] >= 3*np.pi/2) & (sub['elongation_rad'] < 2*np.pi)))
        # Key metric: mean |cosD| (phase coverage quality)
        mean_abs_cosD = float(np.mean(np.abs(sub['cosD'])))
        results[stn] = {
            'n_obs': n,
            'date_range': [date_min, date_max],
            'elongation_mean_rad': elong_mean,
            'elongation_std_rad': elong_std,
            'cosD_mean': cosD_mean,
            'cosD_std': cosD_std,
            'residual_rms_m': rms,
            'phase_coverage': {
                'q1_waxing_crescent': q1,
                'q2_waxing_gibbous': q2,
                'q3_waning_gibbous': q3,
                'q4_waning_crescent': q4
            },
            'mean_abs_cosD': mean_abs_cosD
        }
    return results


def analyze_temporal_epochs(df):
    """Analyze how station composition changes across temporal epochs."""
    epochs = [
        ('pre_1990', df['date_julian_year'] < 1990),
        ('1990_2000', (df['date_julian_year'] >= 1990) & (df['date_julian_year'] < 2000)),
        ('2000_2010', (df['date_julian_year'] >= 2000) & (df['date_julian_year'] < 2010)),
        ('post_2010', df['date_julian_year'] >= 2010)
    ]
    results = {}
    for name, mask in epochs:
        sub = df[mask]
        station_counts = {s: int(np.sum(sub['station'] == s)) for s in sub['station'].unique()}
        results[name] = {
            'n_total': len(sub),
            'station_counts': station_counts,
            'dominant_station': max(station_counts, key=station_counts.get) if station_counts else None,
            'mean_abs_cosD': float(np.mean(np.abs(sub['cosD'])))
        }
    return results


def analyze_covariate_shift(df, split_jd=2454600):
    """Analyze how station mix and elongation distribution differ pre vs post split."""
    pre = df[df['date_julian'] < split_jd]
    post = df[df['date_julian'] >= split_jd]
    
    results = {
        'split_jd': split_jd,
        'pre': {
            'n': len(pre),
            'stations': {s: int(np.sum(pre['station'] == s)) for s in pre['station'].unique()},
            'elongation_mean': float(pre['elongation_rad'].mean()),
            'elongation_std': float(pre['elongation_rad'].std()),
            'cosD_mean': float(pre['cosD'].mean()),
            'cosD_std': float(pre['cosD'].std()),
            'mean_abs_cosD': float(np.mean(np.abs(pre['cosD'])))
        },
        'post': {
            'n': len(post),
            'stations': {s: int(np.sum(post['station'] == s)) for s in post['station'].unique()},
            'elongation_mean': float(post['elongation_rad'].mean()),
            'elongation_std': float(post['elongation_rad'].std()),
            'cosD_mean': float(post['cosD'].mean()),
            'cosD_std': float(post['cosD'].std()),
            'mean_abs_cosD': float(np.mean(np.abs(post['cosD'])))
        }
    }
    
    # KS test for elongation and cosD distributions
    ks_elong = stats.ks_2samp(pre['elongation_rad'].values, post['elongation_rad'].values)
    ks_cosD = stats.ks_2samp(pre['cosD'].values, post['cosD'].values)
    results['ks_test'] = {
        'elongation_statistic': float(ks_elong.statistic),
        'elongation_pvalue': float(ks_elong.pvalue),
        'cosD_statistic': float(ks_cosD.statistic),
        'cosD_pvalue': float(ks_cosD.pvalue)
    }
    return results


def analyze_leave_one_station_out(df):
    """Analyze what happens when each station is held out: covariate shift."""
    stations = df['station'].unique()
    results = {}
    for held_out in stations:
        train = df[df['station'] != held_out]
        test = df[df['station'] == held_out]
        
        # Compare elongation distributions
        ks = stats.ks_2samp(train['elongation_rad'].values, test['elongation_rad'].values)
        
        # Effective leverage: how much does cosD range differ?
        train_cosD_range = float(np.ptp(train['cosD']))
        test_cosD_range = float(np.ptp(test['cosD']))
        
        results[held_out] = {
            'n_train': len(train),
            'n_test': len(test),
            'train_cosD_range': train_cosD_range,
            'test_cosD_range': test_cosD_range,
            'train_mean_abs_cosD': float(np.mean(np.abs(train['cosD']))),
            'test_mean_abs_cosD': float(np.mean(np.abs(test['cosD']))),
            'train_date_range': [float(train['date_julian_year'].min()), float(train['date_julian_year'].max())],
            'test_date_range': [float(test['date_julian_year'].min()), float(test['date_julian_year'].max())],
            'ks_elongation_statistic': float(ks.statistic),
            'ks_elongation_pvalue': float(ks.pvalue)
        }
    return results


def main():
    print_status("Step 052: Station Sample and Distribution Analysis", "TITLE")
    
    df = load_data()
    print_status(f"Loaded {len(df)} observations from {df['station'].nunique()} stations", "INFO")
    
    print_status("Analyzing per-station coverage...", "PROCESS")
    station_coverage = analyze_station_coverage(df)
    
    print_status("Analyzing temporal epochs...", "PROCESS")
    temporal_epochs = analyze_temporal_epochs(df)
    
    print_status("Analyzing covariate shift (pre/post 2008)...", "PROCESS")
    covariate_shift = analyze_covariate_shift(df, split_jd=2454600)
    
    print_status("Analyzing leave-one-station-out covariate shift...", "PROCESS")
    loso_shift = analyze_leave_one_station_out(df)
    
    # Summary: compute per-station cosD-only regression
    print_status("Computing per-station cosD-only regression...", "PROCESS")
    per_station_regression = {}
    for stn in df['station'].unique():
        sub = df[df['station'] == stn]
        y = sub['residual_m'].values
        X = np.column_stack([sub['cosD'].values, np.ones(len(sub))])
        try:
            coeffs = stable_lstsq(X, y)[0]
            with suppress_scipy_array_api_matmul_runtime_warning(), np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                resid = y - X @ coeffs
            if not np.all(np.isfinite(resid)):
                raise RuntimeError(f"Non-finite residuals in per-station regression for {stn}.")
            mse = np.mean(resid**2)
            se = np.sqrt(mse * np.linalg.pinv(X.T @ X, rcond=1e-10, hermitian=True)[0, 0])
            eta = coeffs[0] / ETA_SCALE_FACTOR  # Convert to Nordtvedt eta
            eta_err = se / ETA_SCALE_FACTOR
            snr = abs(eta) / max(eta_err, 1e-20)
            per_station_regression[stn] = {
                'eta': float(eta),
                'eta_err': float(eta_err),
                'snr': float(snr),
                'n': len(sub)
            }
        except Exception as e:
            per_station_regression[stn] = {'error': str(e)}
    
    results = {
        'step_id': 'step_052',
        'status': 'PASS',
        'step': '052_station_distribution_analysis',
        'total_obs': len(df),
        'station_coverage': station_coverage,
        'temporal_epochs': temporal_epochs,
        'covariate_shift_pre_post': covariate_shift,
        'leave_one_station_out_shift': loso_shift,
        'per_station_cosD_regression': per_station_regression
    }
    
    output_path = PROJECT_ROOT / "results/outputs/step_052_station_distribution_analysis.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print_status(f"Results saved to {output_path}", "SUCCESS")
    
    # Print summary
    print_status("\n=== Station Coverage Summary ===", "TITLE")
    for stn, info in station_coverage.items():
        print_status(f"{stn}: n={info['n_obs']}, years={info['date_range'][0]:.1f}-{info['date_range'][1]:.1f}, "
                    f"RMS={info['residual_rms_m']:.3f}m, mean|cosD|={info['mean_abs_cosD']:.3f}", "INFO")
    
    print_status("\n=== Temporal Epoch Composition ===", "TITLE")
    for epoch, info in temporal_epochs.items():
        print_status(f"{epoch}: n={info['n_total']}, dominant={info['dominant_station']}, "
                    f"mean|cosD|={info['mean_abs_cosD']:.3f}", "INFO")
    
    print_status("\n=== Covariate Shift (pre/post 2008) ===", "TITLE")
    print_status(f"Pre: n={covariate_shift['pre']['n']}, mean|cosD|={covariate_shift['pre']['mean_abs_cosD']:.3f}", "INFO")
    print_status(f"Post: n={covariate_shift['post']['n']}, mean|cosD|={covariate_shift['post']['mean_abs_cosD']:.3f}", "INFO")
    print_status(f"KS test elongation: D={covariate_shift['ks_test']['elongation_statistic']:.4f}, "
                f"p={covariate_shift['ks_test']['elongation_pvalue']:.2e}", "INFO")
    print_status(f"KS test cosD: D={covariate_shift['ks_test']['cosD_statistic']:.4f}, "
                f"p={covariate_shift['ks_test']['cosD_pvalue']:.2e}", "INFO")
    
    print_status("\n=== Leave-One-Station-Out Shift ===", "TITLE")
    for stn, info in loso_shift.items():
        print_status(f"Hold out {stn}: train n={info['n_train']}, test n={info['n_test']}, "
                    f"train|cosD|={info['train_mean_abs_cosD']:.3f}, test|cosD|={info['test_mean_abs_cosD']:.3f}, "
                    f"KS_D={info['ks_elongation_statistic']:.4f}", "INFO")
    
    print_status("\n=== Per-Station cosD Regression ===", "TITLE")
    for stn, info in per_station_regression.items():
        if 'error' not in info:
            print_status(f"{stn}: eta={info['eta']:.3e} +/- {info['eta_err']:.3e}, SNR={info['snr']:.2f}sigma, n={info['n']}", "INFO")


if __name__ == '__main__':
    main()
