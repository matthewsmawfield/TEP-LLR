#!/usr/bin/env python3
"""
Step 025: Leverage Temporal Clustering

Identifies high-leverage points from the OLS fit and bins them into temporal
epochs to investigate if the OLS vs Theil-Sen divergence is driven by localized
observational eras or hardware upgrades.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scripts.utils.numerics import stable_lstsq
import pandas as pd
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status


def compute_cooks_distance(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Compute Cook's Distance for each observation."""
    n = len(y)
    p = X.shape[1]
    
    # Fit OLS
    beta = stable_lstsq(X, y)[0]
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        y_pred = X @ beta
    residuals = y - y_pred
    mse = np.sum(residuals**2) / (n - p)
    
    # leverage h_ii = sum_j Q[i,j]^2 for reduced QR of X
    from scripts.utils.numerics import hat_diagonal_from_qr
    leverage = hat_diagonal_from_qr(X)
    
    # Cook's Distance formula
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        cooks_d = (residuals**2 / (p * mse)) * (leverage / ((1 - leverage)**2))
    return cooks_d


def run_leverage_clustering(logger):
    print_status("═══ Starting Leverage Temporal Clustering analysis...", "TITLE")
    print_status("═══ STEP PURPOSE: Identify high-leverage points from OLS fit and bin into temporal epochs to investigate OLS vs Theil-Sen divergence", "INFO")
    print_status("═══ METHOD: Cook's Distance calculation, temporal binning (5-year epochs), station distribution analysis", "INFO")
    print_status("═══ PARAMETERS: Cook's D threshold=4/n, epoch bin size=5 years", "INFO")

    logger.info(">>> Starting Leverage Temporal Clustering analysis...")

    print_status("═══ DATA SUMMARY", "INFO")
    data_file = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not data_file.exists():
        logger.error(f"✗ Data file not found: {data_file}")
        return None

    df = pd.read_csv(data_file)
    n_total = len(df)
    print_status(f"    Dataset: N = {n_total:,} observations", "DATA")
    print_status(f"    Data source: INPOP19a_all_stations_residuals.csv", "DATA")
    
    x_val = np.cos(df['elongation_rad'].values)
    X = np.column_stack([np.ones(n_total), x_val])
    y = df['residual_m'].values

    # Calculate Cook's Distance
    cooks_d = compute_cooks_distance(X, y)
    
    # Define high leverage threshold (standard Cook's distance criterion: 4/n)
    # Source: Cook & Weisberg 1982, standard diagnostic for influential points
    threshold = 4.0 / n_total
    high_lev_mask = cooks_d > threshold
    
    # Time limits
    min_year = int(np.floor(df['date_julian_year'].min()))
    max_year = int(np.ceil(df['date_julian_year'].max()))

    print_status("═══ ANALYSIS TRACE", "INFO")
    print_status(f">>> Calculating Cook's Distance for all observations", "PROCESS")
    print_status(f">>> Identifying high-leverage points (threshold = {threshold:.6e})", "PROCESS")
    print_status(f">>> Binning high-leverage points into 5-year epochs", "PROCESS")
    print_status(f">>> Analyzing station distribution of high-leverage points", "PROCESS")

    # Bin by 5-year epochs
    bins = np.arange(min_year, max_year + 5, 5)
    df['epoch_bin'] = pd.cut(df['date_julian_year'], bins=bins, right=False)
    
    clustering_results = {}
    total_high_lev = int(np.sum(high_lev_mask))
    
    df_high = df[high_lev_mask]
    
    # Count overall distribution vs high-leverage distribution
    for interval in df['epoch_bin'].cat.categories:
        str_interval = f"{int(interval.left)}-{int(interval.right)}"
        
        mask_interval = df['epoch_bin'] == interval
        n_total_in_bin = int(mask_interval.sum())
        n_high_in_bin = int((mask_interval & high_lev_mask).sum())
        
        pct_high = float(n_high_in_bin / total_high_lev) * 100 if total_high_lev > 0 else 0
        expected_pct = float(n_total_in_bin / n_total) * 100
        
        clustering_results[str_interval] = {
            "total_obs": n_total_in_bin,
            "high_lev_obs": n_high_in_bin,
            "percent_of_all_high_lev": pct_high,
            "expected_percent": expected_pct,
            "overrepresentation_factor": float(pct_high / expected_pct) if expected_pct > 0 else 0
        }

    # Station-specific high leverage count
    station_counts = df_high['station'].value_counts()
    station_high = {str(k): int(v) for k, v in station_counts.items()}
    
    # Did Grasse dominate the high-leverage points?
    grasse_domination = station_high.get('Grasse', 0) / total_high_lev > 0.6 if total_high_lev > 0 else False

    results = {
        "step_id": "step_025",
        "status": "PASS",
        "total_observations": n_total,
        "high_leverage_threshold": float(threshold),
        "total_high_leverage_points": total_high_lev,
        "temporal_bins": clustering_results,
        "station_distribution": station_high,
        "conclusions": {
            "grasse_dominated": bool(grasse_domination),
            "meaningful_clustering": any(b["overrepresentation_factor"] > 1.5 for b in clustering_results.values()),
            "explanation": "High-leverage points were binned over 5-year epochs. An overrepresentation factor > 1.5 indicates a specific hardware era or localized anomaly injected large variance into the OLS fit."
        }
    }

    print_status("═══ RESULTS SUMMARY", "INFO")
    print_status(f"    Total High-Leverage Points: {total_high_lev:,}", "CALC")
    print_status(f"    Cook's D threshold: {threshold:.6e}", "CALC")
    print_status(f"    Clustered Anomalies Detected: {results['conclusions']['meaningful_clustering']}", "CALC")
    print_status(f"    Grasse Dominated: {results['conclusions']['grasse_dominated']}", "CALC")

    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    High-leverage points binned into 5-year epochs", "INFO")
    print_status(f"    Overrepresentation factor > 1.5 indicates specific hardware era or localized anomaly", "INFO")
    print_status(f"    Station distribution analysis tests for instrument-specific effects", "INFO")

    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_025_leverage_temporal_clustering.json", "INFO")
    print_status(f"    Cook's D threshold: 4/n", "INFO")
    print_status(f"    Epoch bin size: 5 years", "INFO")
    print_status(f"    Data source: INPOP19a_all_stations_residuals.csv", "INFO")

    logger.info(f"    Total High-Leverage Points: {total_high_lev}")
    logger.info(f"    Clustered Anomalies Detected? {results['conclusions']['meaningful_clustering']}")
    logger.info(f"✓   Leverage Temporal Clustering Complete.")

    return results

def main():
    # Setup TEPLogger
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_025", str(log_dir / "step_025_leverage_temporal_clustering.log"))
    set_step_logger(logger)
    
    results = run_leverage_clustering(logger)
    
    # Save output to JSON
    output_path = PROJECT_ROOT / "results" / "outputs" / "step_025_leverage_temporal_clustering.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)
    output_rel = output_path.relative_to(PROJECT_ROOT) if output_path.is_relative_to(PROJECT_ROOT) else output_path
    logger.info(f"    Saved results to {output_rel}")

if __name__ == "__main__":
    main()