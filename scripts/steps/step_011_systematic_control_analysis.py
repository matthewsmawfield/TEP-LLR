#!/usr/bin/env python3
"""
Step 011: Systematic Control Analysis for TEP-LLR
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import pandas as pd
import numpy as np
from scipy import stats
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

# Add the project root to the Python path

def compute_partial_correlation(x, y, z, verbose=False):
    # Remove means
    x_resid = x - np.mean(x)
    y_resid = y - np.mean(y)
    z_resid = z - np.mean(z)

    # Regress x on z
    beta_xz = np.sum(x_resid * z_resid) / np.sum(z_resid**2)
    x_pred = beta_xz * z_resid
    x_residuals = x_resid - x_pred

    # Regress y on z
    beta_yz = np.sum(y_resid * z_resid) / np.sum(z_resid**2)
    y_pred = beta_yz * z_resid
    y_residuals = y_resid - y_pred

    # Partial correlation is correlation of residuals
    r_partial, p_partial = stats.pearsonr(x_residuals, y_residuals)

    if verbose:
        print_status("  [CALC] Partial correlation computation:", "CALC")
        print_status(f"  [CALC]    β_xz (x on z): {beta_xz:.6e}", "CALC")
        print_status(f"  [CALC]    β_yz (y on z): {beta_yz:.6e}", "CALC")
        print_status(
            f"  [CALC]    Var(x_resid):  {np.var(x_residuals):.6e}", "CALC")
        print_status(
            f"  [CALC]    Var(y_resid):  {np.var(y_residuals):.6e}", "CALC")
        print_status(f"  [CALC]    r_partial:     {r_partial:.6e}", "CALC")
        print_status(f"  [CALC]    p_partial:     {p_partial:.2e}", "CALC")

    return r_partial, p_partial, x_residuals, y_residuals

def run_control_analysis(df, verbose=False):
    print_status("="*60, "INFO")
    print_status("SYSTEMATIC CONTROL ANALYSIS - DETAILED TRACE", "TITLE")
    print_status("="*60, "INFO")

    # Extract data
    residuals = df['residual_m'].values
    cos_elong = np.cos(df['elongation_rad'].values)
    jd = df['date_julian'].values
    n = len(residuals)

    print_status(f"[DATA] Dataset: N={n:,} observations", "INFO")
    print_status(
        f"[DATA] Residual range: [{np.min(residuals):.6f}, {np.max(residuals):.6f}] m", "INFO")
    print_status(
        f"[DATA] cos(elongation) range: [{np.min(cos_elong):.6f}, {np.max(cos_elong):.6f}]", "INFO")
    print_status(
        f"[DATA] JD range: [{np.min(jd):.4f}, {np.max(jd):.4f}]", "INFO")

    # Original correlation (no controls)
    r_orig, p_orig = stats.pearsonr(residuals, cos_elong)
    snr_orig = abs(r_orig)/np.sqrt((1-r_orig**2)/(n-2))

    print_status("", "INFO")
    print_status("TEST 1: UNCONTROLLED CORRELATION (BASELINE)", "PROCESS")
    print_status(f"  [CALC] r_original = {r_orig:.6e}", "CALC")
    print_status(f"  [CALC] p_original = {p_orig:.2e}", "CALC")
    print_status(f"  [CALC] Significance: {snr_orig:.2f}σ", "CALC")

    # Control 1: Linear time trend
    jd_norm = (jd - np.mean(jd)) / np.std(jd)
    r_time, p_time, _, _ = compute_partial_correlation(
        residuals, cos_elong, jd_norm, verbose=verbose)

    print_status("", "INFO")
    print_status("TEST 2: CONTROLLING FOR LINEAR TIME TREND", "PROCESS")
    print_status(f"  [CALC] Time trend coefficient: {np.corrcoef(residuals, jd_norm)[0,1]:.6e}", "CALC")
    print_status(f"  [CALC] r_partial (time controlled): {r_time:.6e}", "CALC")
    print_status(f"  [CALC] p_partial: {p_time:.2e}", "CALC")
    print_status(f"  [CALC] Signal attenuation: {(1 - abs(r_time)/abs(r_orig))*100:.1f}%", "CALC")

    # Control 2: Quadratic time trend
    jd_norm_sq = jd_norm**2
    r_quad, p_quad, _, _ = compute_partial_correlation(
        residuals, cos_elong, jd_norm_sq, verbose=verbose)

    print_status("", "INFO")
    print_status("TEST 3: CONTROLLING FOR QUADRATIC TIME TREND", "PROCESS")
    print_status(f"  [CALC] r_partial (quad controlled): {r_quad:.6e}", "CALC")
    print_status(f"  [CALC] p_partial: {p_quad:.2e}", "CALC")

    # Control 3: Seasonal effects (annual cycle)
    if 'date_julian_year' in df.columns:
        t_year = df['date_julian_year'].values
        year_frac = t_year - np.floor(t_year)
        sin_year = np.sin(2 * np.pi * year_frac)
        cos_year = np.cos(2 * np.pi * year_frac)
        
        # Control 4: Lunar Sidereal Cycle (27.32 days)
        # Control for lunar orientation artifacts
        sidereal_period = 27.321661
        sidereal_frac = (jd / sidereal_period) % 1.0
        sin_sidereal = np.sin(2 * np.pi * sidereal_frac)
        cos_sidereal = np.cos(2 * np.pi * sidereal_frac)
        
        # Control 5: Lunar Anomalistic Cycle (27.55 days)
        # Control for distance-dependent artifacts (perigee/apogee)
        anomalistic_period = 27.554551
        anomalistic_frac = (jd / anomalistic_period) % 1.0
        sin_anomalistic = np.sin(2 * np.pi * anomalistic_frac)
        cos_anomalistic = np.cos(2 * np.pi * anomalistic_frac)

        # Combined control: all systematics simultaneously including lunar cycles
        X = np.column_stack([
            np.ones(n), 
            jd_norm, 
            jd_norm**2,
            sin_year, 
            cos_year,
            sin_sidereal,
            cos_sidereal,
            sin_anomalistic,
            cos_anomalistic
        ])

        # Regress residuals on controls
        beta_res = np.linalg.lstsq(X, residuals, rcond=None)[0]
        res_residuals = residuals - X @ beta_res

        # Regress cos_elong on controls
        beta_cos = np.linalg.lstsq(X, cos_elong, rcond=None)[0]
        cos_residuals = cos_elong - X @ beta_cos

        # Partial correlation after all controls
        r_combined, p_combined = stats.pearsonr(res_residuals, cos_residuals)
        
        print_status("", "INFO")
        print_status("TEST 4: CONTROLLING FOR SEASONAL AND LUNAR CYCLES", "PROCESS")
        print_status(f"  [CALC] Controlled for: Annual, Sidereal (27.32d), Anomalistic (27.55d)", "CALC")
        print_status(f"  [CALC] r_combined (all controlled): {r_combined:.6e}", "CALC")
        print_status(f"  [CALC] p_combined: {p_combined:.2e}", "CALC")
    else:
        r_combined, p_combined = r_orig, p_orig

    if verbose:
        print_status("", "INFO")
        print_status("TEST 5: COMBINED CONTROL (ALL SYSTEMATICS)", "PROCESS")
        print_status(
            "  [CALC] Multiple regression coefficients for residuals:", "CALC")
        print_status(f"  [CALC]    β_0 (intercept): {beta_res[0]:.6e}", "CALC")
        print_status(
            f"  [CALC]    β_1 (linear):      {beta_res[1]:.6e}", "CALC")
        print_status(
            f"  [CALC]    β_2 (quadratic):   {beta_res[2]:.6e}", "CALC")
        print_status(
            f"  [CALC]    Var(res_residuals): {np.var(res_residuals):.6e}", "CALC")
        print_status(
            f"  [CALC]    Var(cos_residuals): {np.var(cos_residuals):.6e}", "CALC")
        print_status(f"  [CALC] r_combined: {r_combined:.6e}", "CALC")
        print_status(f"  [CALC] p_combined: {p_combined:.2e}", "CALC")
        print_status(
            f"  [CALC] Total attenuation: {(1 - abs(r_combined)/abs(r_orig))*100:.1f}%", "CALC")
        print_status("", "INFO")
        print_status("="*60, "INFO")
        final_ok = bool(p_combined < 0.05)
        final_msg = "FINAL RESULT: Signal persists after all systematic controls" if final_ok else "FINAL RESULT: Signal not robust after all systematic controls"
        print_status(final_msg, "SUCCESS" if final_ok else "WARNING")
        print_status(f"  r_original:   {r_orig:.6e} (p={p_orig:.2e})", "CALC")
        print_status(
            f"  r_combined:   {r_combined:.6e} (p={p_combined:.2e})", "CALC")
        print_status("="*60, "INFO")

    return {
        "r_original": float(r_orig),
        "p_original": float(p_orig),
        "r_partial_linear_time": float(r_time),
        "p_partial_linear_time": float(p_time),
        "r_partial_quadratic_time": float(r_quad),
        "p_partial_quadratic_time": float(p_quad),
        "r_partial_seasonal": float(r_combined) if 'r_combined' in locals() else None,
        "p_partial_seasonal": float(p_combined) if 'p_combined' in locals() else None,
        "r_partial_all_controls": float(r_combined),
        "p_partial_all_controls": float(p_combined),
        "signal_persists": bool(p_combined < 0.05),
        "attenuation_percent": float((1 - abs(r_combined)/abs(r_orig))*100)
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 011: Systematic Control Analysis")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_011", str(
        log_dir / "step_011_systematic_control_analysis.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Starting Systematic Control Analysis...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    summary = run_control_analysis(df)

    results = {
        "step_id": "step_011",
        "control_results": summary,
        "status": "PASS" if summary["signal_persists"] else "WARNING"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_011_systematic_control_analysis")
    print_status("Control Analysis Complete.", "SUCCESS")