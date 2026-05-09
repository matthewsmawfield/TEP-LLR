#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 052: Modern-Era Power Design Analysis

Computes the power required to detect a signal of amplitude η = -3 × 10⁻⁴ using:
N_required = (13|η|σ_cosD / (z_ασ_r))²

Reports power table for subsets: INPOP full, C-SPAD, APO only, Grasse only, DE430
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import detect_outliers_sigma

TEP_CONFIG = get_config()

def calculate_power_stats(df_subset, target_eta=-3e-4, z_alpha=3.0, apply_outlier_cleaning=False):
    """
    Calculate power statistics for a data subset.
    
    Args:
        df_subset: DataFrame with residual data
        target_eta: Target signal amplitude for power calculation
        z_alpha: Critical value for significance (3 for 3σ)
        apply_outlier_cleaning: Whether to apply 6σ MAD outlier cleaning (default False for consistency with existing behavior)
    
    Returns:
        dict with power statistics
    """
    df_work = df_subset.copy()  # Keep this copy since we modify it
    outlier_info = {
        'outlier_cleaning_applied': apply_outlier_cleaning,
        'n_outliers_removed': 0,
    }

    if apply_outlier_cleaning:
        outlier_mask = detect_outliers_sigma(df_work['residual_m'].values, sigma_threshold=6.0)
        n_outliers = int(np.sum(outlier_mask))
        df_work = df_work[~outlier_mask]  # PERFORMANCE FIX: Removed unnecessary .copy()
        outlier_info['n_outliers_removed'] = n_outliers
        outlier_info['n_cleaned'] = len(df_work)
    
    N = len(df_work)
    sigma_r = df_work['residual_m'].std()
    cos_D = np.cos(df_work['elongation_rad'].values)
    sigma_cosD = cos_D.std()
    
    # Calculate expected SNR
    expected_snr = (13.0 * abs(target_eta) * sigma_cosD * np.sqrt(N)) / (z_alpha * sigma_r)
    
    # Calculate N_required for 3σ detection
    # N_required = (z_alpha * sigma_r / (13 * |eta| * sigma_cosD))^2
    N_required = ((z_alpha * sigma_r) / (13.0 * abs(target_eta) * sigma_cosD))**2
    
    # Estimate η from simple regression
    # polyfit returns amplitude in meters, need to convert to eta by dividing by ETA_SCALE_FACTOR
    eta_est = np.polyfit(cos_D, df_work['residual_m'].values, 1)[0] / ETA_SCALE_FACTOR
    
    # Phase coverage
    phase_range = df_work['elongation_rad'].max() - df_work['elongation_rad'].min()
    phase_coverage = f"{phase_range/np.pi:.1f}π"
    
    result = {
        'N': int(N),
        'RMS_cm': float(sigma_r * 100),  # Convert to cm
        'phase_coverage': phase_coverage,
        'expected_snr': float(expected_snr),
        'N_required': float(N_required),
        'observed_eta': float(eta_est)
    }
    
    # Add outlier cleaning info if applied
    if apply_outlier_cleaning:
        result['outlier_cleaning'] = outlier_info
    
    return result

def run_power_design_analysis(verbose=False):
    """Run power design analysis for all data subsets."""
    input_path_inpop = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    input_path_de430 = PROJECT_ROOT / "data" / "processed" / "DE430_all_residuals.csv"
    
    if not input_path_inpop.exists():
        print_status(f"CRITICAL DATA FAILURE: {input_path_inpop} not found. Cannot proceed.", "ERROR")
        return None
    
    df_inpop = pd.read_csv(input_path_inpop)
    
    print_status("Running power design analysis...", "INFO")
    
    # Add year column if not present
    if 'year' not in df_inpop.columns:
        df_inpop['year'] = df_inpop['date_julian_year'].values
    
    target_eta = -3e-4
    
    # Analyze subsets
    subsets = {}
    
    # INPOP full
    subsets['INPOP full'] = calculate_power_stats(df_inpop, target_eta=target_eta)
    
    # C-SPAD (2009-2019)
    df_cspad = df_inpop[(df_inpop['year'] >= 2009) & (df_inpop['year'] <= 2019)]
    subsets['C-SPAD (2009-2019)'] = calculate_power_stats(df_cspad, target_eta=target_eta)
    
    # APO only
    df_apo = df_inpop[df_inpop['station'] == 'APO']
    subsets['APO only'] = calculate_power_stats(df_apo, target_eta=target_eta)
    
    # Grasse only
    df_grasse = df_inpop[df_inpop['station'] == 'Grasse']
    subsets['Grasse only'] = calculate_power_stats(df_grasse, target_eta=target_eta)
    
    # DE430 (if available) - CRITICAL FIX: Apply outlier cleaning to match step_005 and manuscript
    if input_path_de430.exists():
        df_de430 = pd.read_csv(input_path_de430)
        print_status("Applying 6σ MAD outlier cleaning to DE430 (CRITICAL FIX)", "INFO")
        subsets['DE430'] = calculate_power_stats(df_de430, target_eta=target_eta, apply_outlier_cleaning=True)
    else:
        print_status("DE430 data not available, skipping DE430 subset.", "WARNING")
    
    # Print results
    for subset_name, stats in subsets.items():
        print_status(f"  {subset_name}:", "INFO")
        print_status(f"    N = {stats['N']:,}, RMS = {stats['RMS_cm']:.1f} cm", "INFO")
        print_status(f"    Expected SNR = {stats['expected_snr']:.2f}σ, Observed η = {stats['observed_eta']:.8e}", "INFO")
    
    results = {
        'target_eta': float(target_eta),
        'z_alpha': 3.0,
        'subsets': subsets
    }
    
    return results

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_052", str(log_dir / "step_052_power_design_analysis.log"))
    set_step_logger(logger)
    
    results = run_power_design_analysis(verbose=True)
    
    if results:
        logger.save_step_results(results, PROJECT_ROOT, "step_052_power_design_analysis")
        print_status("Power Design Analysis Complete.", "SUCCESS")
    else:
        print_status("Power Design Analysis Failed.", "ERROR")
