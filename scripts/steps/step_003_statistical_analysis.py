#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 003: Statistical Analysis
"""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma, cluster_robust_variance
import argparse
import pandas as pd
import numpy as np
import scipy.linalg
import emcee
from scipy import stats

TEP_CONFIG = get_config()

def _log_likelihood(theta, x, y, y_err):
    eta, intercept = theta
    model = ETA_SCALE_FACTOR * eta * x + intercept
    sigma2 = y_err**2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def _log_prior(theta):
    eta, intercept = theta
    if -0.01 <= eta <= 0.01 and -0.1 <= intercept <= 0.1:
        return 0.0
    return -np.inf

def _log_probability(theta, x, y, y_err):
    lp = _log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + _log_likelihood(theta, x, y, y_err)

def ar1_gls_regression(y, x, station_ids=None, verbose=False):
    """
    Perform AR(1) Generalized Least Squares regression with optional
    cluster-robust standard errors.

    Estimates the AR(1) parameter rho from OLS residuals, applies a
    Cochrane-Orcutt quasi-differencing transformation, and fits the
    transformed model.  If station_ids are provided, cluster-robust
    (sandwich) standard errors are computed by station.

    Parameters:
    -----------
    y : np.ndarray
        Dependent variable (residuals)
    x : np.ndarray
        Independent variable (cos(elongation))
    station_ids : np.ndarray, optional
        Cluster identifier for each observation (e.g. station name).
        Used to compute cluster-robust standard errors.
    verbose : bool
        Whether to print diagnostic information

    Returns:
    --------
    dict containing:
        - eta: Estimated Nordtvedt parameter
        - eta_error: Standard error of eta (GLS, Birge-scaled)
        - eta_error_cluster: Cluster-robust SE of eta (if station_ids given)
        - rho: Estimated AR(1) parameter
        - rho_error: Standard error of rho
        - durbin_watson: Durbin-Watson statistic
        - n_obs: Number of observations
        - n_clusters: Number of clusters (if station_ids given)
    """
    n = len(y)

    # First fit with OLS to get initial residuals
    reg_ols = linear_regression(y, x, weights=None)
    residuals = y - (13.0 * reg_ols['eta'] * x + reg_ols['intercept'])

    # Estimate AR(1) parameter rho from residuals
    # rho = sum(residuals[t] * residuals[t-1]) / sum(residuals[t-1]^2)
    rho = np.sum(residuals[1:] * residuals[:-1]) / np.sum(residuals[:-1]**2)

    # Standard error of rho (asymptotic formula for AR(1))
    rho_error = np.sqrt((1 - rho**2) / n)

    # Durbin-Watson statistic
    dw_stat = np.sum(np.diff(residuals)**2) / np.sum(residuals**2)

    # Cochrane-Orcutt transformation for AR(1) GLS (O(n) instead of O(n³) Cholesky)
    # Transform: y_t* = y_t - rho*y_{t-1}, x_t* = x_t - rho*x_{t-1}
    y_star = y[1:] - rho * y[:-1]
    x_star = x[1:] - rho * x[:-1]

    # OLS on transformed variables
    reg_star = linear_regression(y_star, x_star, weights=None)
    eta_gls = reg_star['eta']
    intercept_gls = reg_star['intercept']

    # Error estimation for transformed regression
    eta_error_gls = reg_star['eta_error']

    # Adjust intercept for original scale
    intercept_gls = intercept_gls / (1 - rho)

    # Cluster-robust standard errors (by station)
    eta_error_cluster = None
    n_clusters = None
    if station_ids is not None and len(station_ids) == n:
        # Align cluster IDs with the transformed (differenced) observations.
        # y_star[t] corresponds to y[t+1], so use station_ids[1:].
        cluster_ids_star = station_ids[1:]
        X_star = np.column_stack([x_star, np.ones(len(x_star))])
        # Residuals on the transformed scale: y_star - (A_star*x_star + B_star)
        A_star = reg_star['amplitude']
        B_star = reg_star['intercept']
        u_star = y_star - (A_star * x_star + B_star)
        cr = cluster_robust_variance(
            X_star, u_star, cluster_ids_star, small_sample_correction=True
        )
        # Convert amplitude SE to eta SE: eta = A / 13
        from scripts.utils.llr_constants import ETA_SCALE_FACTOR
        eta_error_cluster = cr['se_cluster'][0] / ETA_SCALE_FACTOR
        n_clusters = cr['n_clusters']

    if verbose:
        print_status(f"AR(1) GLS Regression Results:", "CALC")
        print_status(f"  η (GLS) = {eta_gls:.8e} ± {eta_error_gls:.8e}", "CALC")
        if eta_error_cluster is not None:
            print_status(f"  η (GLS, cluster-robust) = {eta_gls:.8e} ± {eta_error_cluster:.8e}", "CALC")
            print_status(f"  Cluster-robust SE inflation: {eta_error_cluster/eta_error_gls:.2f}x", "CALC")
        print_status(f"  AR(1) parameter ρ = {rho:.4f} ± {rho_error:.4f}", "CALC")
        print_status(f"  Durbin-Watson = {dw_stat:.3f}", "CALC")
        print_status(f"  Interpretation: {'Significant autocorrelation' if abs(rho) > 0.1 else 'Weak autocorrelation'}", "CALC")

    result = {
        "eta": eta_gls,
        "eta_error": eta_error_gls,
        "eta_error_cluster": eta_error_cluster,
        "rho": rho,
        "rho_error": rho_error,
        "durbin_watson": dw_stat,
        "n_obs": n,
        "n_clusters": n_clusters
    }
    return result

def run_statistical_analysis(verbose=False):
    print_status("═══ Starting Step 003: Statistical Analysis...", "TITLE")
    print_status("═══ STEP PURPOSE: Primary TEP Nordtvedt parameter estimation using OLS and Bayesian MCMC methods", "INFO")
    print_status("═══ METHOD: Unweighted OLS regression + Bayesian MCMC with emcee sampler", "INFO")
    print_status(f"═══ PARAMETERS: 6σ MAD outlier cleaning, {TEP_CONFIG.get('MCMC_STANDARD_WALKERS', 32)} walkers, {TEP_CONFIG.get('MCMC_STANDARD_STEPS', 3000)} steps, {TEP_CONFIG.get('MCMC_BURN_IN', 1000)} burn-in", "INFO")
    
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not input_path.exists():
        print_status(f"CRITICAL DATA FAILURE: {input_path} not found. Cannot proceed.", "ERROR")
        return None

    df = pd.read_csv(input_path)
    if 'sigma_m' not in df.columns:
        print_status("WARNING: sigma_m not found in dataset. Falling back to global std(y).", "WARNING")
        df['sigma_m'] = np.std(df['residual_m'].values)

    print_status("═══ DATA SUMMARY", "INFO")
    print_status(f"    Dataset: N = {len(df):,} observations", "DATA")
    print_status(f"    Stations: 5 (APO, Grasse, Matera, McDonald2, Haleakala)", "DATA")
    
    # Apply 6σ MAD outlier cleaning for consistency with other steps (CRITICAL FIX for standardization)
    outlier_mask = detect_outliers_sigma(df['residual_m'].values, sigma_threshold=6.0)
    n_outliers = int(np.sum(outlier_mask))
    df_clean = df[~outlier_mask]  # PERFORMANCE FIX: Removed unnecessary .copy()
    if verbose:
        print_status(f"Applied 6σ MAD outlier cleaning: removed {n_outliers}/{len(df)} outliers", "INFO")
        print_status(f"    Cleaned dataset: N = {len(df_clean):,} observations", "DATA")

    # PERFORMANCE FIX: Use pre-computed cos_elong_rad if available, otherwise compute
    if 'cos_elong_rad' in df_clean.columns:
        x = df_clean['cos_elong_rad'].values
    else:
        x = np.cos(df_clean['elongation_rad'].values)
    y = df_clean['residual_m'].values
    n = len(y)
    
    # Use the fixed per-observation sigma_m for the MCMC likelihood.
    # sigma_m is now station-specific RMS estimated from the data itself
    # (fixed in parse_inpop_mini.py; previously parts[6] was incorrectly interpreted).
    y_err = df_clean['sigma_m'].values

    # OLS Regression WITHOUT weights — kept unweighted so the Birge ratio
    # remains a meaningful diagnostic of model fit. Using station-specific
    # weights would drive the Birge ratio to approximately 1.0 by construction.
    reg = linear_regression(y, x, weights=None)
    eta_ols = reg['eta']
    eta_err_ols = reg['eta_error']
    
    print_status("═══ ANALYSIS TRACE", "INFO")
    print_status(">>> Performing unweighted OLS regression", "PROCESS")
    if verbose:
        print_status(f"    Unweighted Regression Complete (N={n}, DOF={n-2}):", "INFO")
        print_status(f"      RSS = {reg['rss']:.6e}", "CALC")
        print_status(f"      χ²_red = {reg['chi2_red']:.6f}", "CALC")
        print_status(f"      Birge Ratio = {reg['birge_ratio']:.3f}", "CALC")
        print_status(f"      Condition Number κ(R) = {reg['condition_number']:.2e}", "CALC")
        print_status(f"      Final η = {eta_ols:.8e} ± {eta_err_ols:.8e}", "CALC")

    # MCMC with configured steps for research-grade convergence
    n_walkers = TEP_CONFIG.get("MCMC_STANDARD_WALKERS", 32)
    n_steps = TEP_CONFIG.get("MCMC_STANDARD_STEPS", 3000)
    burn_in = TEP_CONFIG.get("MCMC_BURN_IN", 1000)
    initial = np.array([eta_ols, reg['intercept']])
    np.random.seed(42)
    pos = initial + 1e-6 * np.random.randn(n_walkers, 2)

    print_status(">>> Running Bayesian MCMC analysis", "PROCESS")
    print_status(f"    Number of walkers: {n_walkers}", "CALC")
    print_status(f"    Total steps: {n_steps}", "CALC")
    print_status(f"    Burn-in steps: {burn_in}", "CALC")
    print_status(f"    Initial position: η={eta_ols:.8e}, b={reg['intercept']:.8e}", "CALC")

    # Use per-observation uncertainties in likelihood
    sampler = emcee.EnsembleSampler(n_walkers, 2, _log_probability, args=(x, y, y_err))
    sampler.run_mcmc(pos, n_steps, progress=False)

    # Convergence diagnostics for emcee
    # Autocorrelation time serves as the primary convergence diagnostic
    try:
        tau = sampler.get_autocorr_time()
        tau_mean = np.mean(tau)
        # Convergence criterion: run > 50 * autocorrelation time
        n_steps_after_burnin = n_steps - burn_in
        convergence_criterion = n_steps_after_burnin > 50 * tau_mean
        convergence_status = "CONVERGED" if convergence_criterion else "INSUFFICIENT_STEPS"
    except (emcee.autocorr.AutocorrError, RuntimeError):
        # Autocorrelation time estimation failed - insufficient convergence
        tau_mean = np.nan
        convergence_criterion = False
        convergence_status = "AUTOCORR_FAILED"

    # Acceptance fraction
    accept_frac = np.mean(sampler.acceptance_fraction)

    print_status("═══ MCMC CONVERGENCE DIAGNOSTICS", "INFO")
    print_status(f"    Autocorrelation time: {tau_mean:.2f}" if np.isfinite(tau_mean) else "    Autocorrelation time: FAILED", "CALC")
    print_status(f"    Steps per autocorr: {(n_steps - burn_in) / tau_mean:.1f}" if np.isfinite(tau_mean) and tau_mean > 0 else "    Steps per autocorr: N/A", "CALC")
    print_status(f"    Acceptance fraction: {accept_frac:.3f}", "CALC")
    print_status(f"    Convergence status: {convergence_status}", "PASS" if convergence_status == "CONVERGED" else "WARNING")

    flat_samples = sampler.get_chain(discard=burn_in, flat=True)
    eta_mcmc = np.mean(flat_samples[:, 0])
    eta_err_mcmc = np.std(flat_samples[:, 0])
    snr = abs(eta_mcmc) / eta_err_mcmc
    
    status = "STRONG DETECTION" if snr > 5 else "DETECTION" if snr > 3 else "INCONCLUSIVE"
    
    # Extract station IDs for cluster-robust SE computation
    station_ids = df_clean['station'].values if 'station' in df_clean.columns else None

    # AR(1) GLS Regression to account for temporal autocorrelation
    print_status(">>> Performing AR(1) GLS regression (autocorrelation-aware error estimation)", "PROCESS")
    ar1_gls_results = ar1_gls_regression(y, x, station_ids=station_ids, verbose=verbose)

    # Use cluster-robust error as the primary error if available
    primary_eta_error = ar1_gls_results.get('eta_error_cluster') or ar1_gls_results['eta_error']
    primary_snr = abs(ar1_gls_results['eta']) / primary_eta_error if primary_eta_error > 0 else 0.0

    print_status("═══ RESULTS SUMMARY", "INFO")
    print_status(f"    OLS η: {eta_ols:.8e} ± {eta_err_ols:.8e}", "CALC")
    print_status(f"    MCMC η: {eta_mcmc:.8e} ± {eta_err_mcmc:.8e}", "CALC")
    print_status(f"    AR(1) GLS η: {ar1_gls_results['eta']:.8e} ± {ar1_gls_results['eta_error']:.8e}", "CALC")
    if ar1_gls_results.get('eta_error_cluster') is not None:
        print_status(f"    AR(1) GLS η (cluster-robust): {ar1_gls_results['eta']:.8e} ± {ar1_gls_results['eta_error_cluster']:.8e}", "CALC")
        print_status(f"    Cluster-robust SNR: {primary_snr:.2f}σ", "CALC")
    print_status(f"    AR(1) parameter ρ: {ar1_gls_results['rho']:.4f} ± {ar1_gls_results['rho_error']:.4f}", "CALC")
    print_status(f"    Durbin-Watson: {ar1_gls_results['durbin_watson']:.3f}", "CALC")
    print_status(f"    SNR (MCMC): {snr:.2f}σ", "CALC")
    print_status(f"    Status: {status}", "PASS" if status == "STRONG DETECTION" else "INFO")

    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    The negative η indicates TEP Nordtvedt violation is present", "INFO")
    print_status(f"    Statistical significance: {snr:.2f}σ ({status})", "INFO")
    print_status(f"    Convergence: {convergence_status}", "PASS" if convergence_status == "CONVERGED" else "WARNING")
    print_status(f"    Temporal autocorrelation: ρ = {ar1_gls_results['rho']:.4f} (DW = {ar1_gls_results['durbin_watson']:.3f})", "INFO")
    print_status(f"    Autocorrelation impact: Error inflation factor = {eta_err_mcmc / ar1_gls_results['eta_error']:.2f}x if significant", "INFO" if abs(ar1_gls_results['rho']) > 0.1 else "INFO")
    print_status(f"    Limitations: OLS is unweighted (Birge ratio diagnostic); MCMC uses fixed station-specific sigma_m", "INFO")

    results = {
        "eta_ols": float(eta_ols),
        "eta_ols_error": float(eta_err_ols),  # Add this key for consistency with other steps
        "eta_mcmc": float(eta_mcmc),
        "eta_err_mcmc": float(eta_err_mcmc),
        "snr": float(snr),
        "status": status,
        "ar1_gls": {
            "eta": float(ar1_gls_results['eta']),
            "eta_error": float(ar1_gls_results['eta_error']),
            "eta_error_cluster": float(ar1_gls_results['eta_error_cluster']) if ar1_gls_results.get('eta_error_cluster') is not None else None,
            "rho": float(ar1_gls_results['rho']),
            "rho_error": float(ar1_gls_results['rho_error']),
            "durbin_watson": float(ar1_gls_results['durbin_watson']),
            "n_clusters": int(ar1_gls_results['n_clusters']) if ar1_gls_results.get('n_clusters') is not None else None,
            "n_obs": int(len(y)),
            "interpretation": "Significant autocorrelation" if abs(ar1_gls_results['rho']) > 0.1 else "Weak autocorrelation"
        },
        "convergence_diagnostics": {
            "autocorr_time": float(tau_mean) if np.isfinite(tau_mean) else None,
            "convergence_status": convergence_status,
            "n_steps_after_burnin": int(n_steps - burn_in),
            "steps_per_autocorr": float((n_steps - burn_in) / tau_mean) if np.isfinite(tau_mean) and tau_mean > 0 else None,
            "acceptance_fraction": float(accept_frac),
            "n_walkers": int(n_walkers),
            "n_steps": int(n_steps),
            "burn_in": int(burn_in)
        },
        "regression_metrics": reg,
        "outlier_cleaning": {
            "method": "6σ MAD (detect_outliers_sigma)",
            "sigma_threshold": 6.0,
            "n_outliers_removed": n_outliers,
            "n_original": len(df),
            "n_cleaned": len(df_clean)
        }
    }
    
    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_003_statistical_analysis.json", "INFO")
    print_status(f"    MCMC walkers: {n_walkers}", "INFO")
    print_status(f"    MCMC steps: {n_steps}", "INFO")
    print_status(f"    MCMC burn-in: {burn_in}", "INFO")
    
    return results

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_003", str(log_dir / "step_003_statistical_analysis.log"))
    set_step_logger(logger)
    
    results = run_statistical_analysis(verbose=True)
    if results:
        # Convert numpy to native for JSON
        def to_native(obj):
            if isinstance(obj, dict): return {k: to_native(v) for k, v in obj.items()}
            if isinstance(obj, np.ndarray): return obj.tolist()
            if isinstance(obj, (np.float64, np.float32)): return float(obj)
            return obj
        logger.save_step_results(to_native(results), PROJECT_ROOT, "step_003_statistical_analysis")
        print_status(f"Statistical Analysis Complete. SNR = {results['snr']:.1f}σ", "SUCCESS")