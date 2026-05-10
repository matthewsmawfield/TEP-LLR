#!/usr/bin/env python3
"""
Step 016: Bayesian MCMC Analysis for TEP-LLR
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import numpy as np
import pandas as pd
import emcee
from scipy import stats

from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.statistical_utils import detect_outliers_sigma

TEP_CONFIG = get_config()


def gelman_rubin(chains):
    """Compute Gelman-Rubin statistic for convergence diagnostics."""
    n_walkers, n_steps, n_params = chains.shape

    # Reshape to (n_walkers, n_steps, n_params)
    # Compute within-chain variance
    W = np.mean(np.var(chains, axis=1, ddof=1), axis=0)

    # Compute between-chain variance
    chain_means = np.mean(chains, axis=1)
    B = n_steps * np.var(chain_means, axis=0, ddof=1)

    # Estimated variance
    V_hat = ((n_steps - 1) / n_steps) * W + (1 / n_steps) * B

    # Gelman-Rubin statistic
    R_hat = np.sqrt(V_hat / W)

    return R_hat


# Alias for compatibility
compute_gelman_rubin = gelman_rubin


def log_probability(theta, x, y, y_err):
    """Log probability function for MCMC with uniform prior bounds."""
    eta, intercept = theta

    # CRITICAL FIX: enforce prior bounds [-0.01, 0.01] for eta and [-0.1, 0.1] for intercept.
    # These bounds were stated in verbose output but never actually applied to the sampler,
    # meaning the Savage-Dickey Bayes factor was computed under a prior the walkers did not respect.
    if not (-0.01 <= eta <= 0.01 and -0.1 <= intercept <= 0.1):
        return -np.inf

    model = eta * ETA_SCALE_FACTOR * x + intercept
    residuals = y - model
    chi2 = np.sum((residuals / y_err) ** 2)
    return -0.5 * chi2

def run_bayesian_analysis(verbose=False):
    if verbose:
        print_status("="*60, "INFO")
        print_status("BAYESIAN MCMC ANALYSIS - DETAILED TRACE", "TITLE")
        print_status("="*60, "INFO")

    input_path = PROJECT_ROOT / "data" / "processed" / \
        "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)

    # Apply 6σ MAD outlier cleaning for consistency with other steps (CRITICAL FIX for standardization)
    outlier_mask = detect_outliers_sigma(df['residual_m'].values, sigma_threshold=6.0)
    n_outliers = int(np.sum(outlier_mask))
    df_clean = df[~outlier_mask]  # PERFORMANCE FIX: Removed unnecessary .copy()
    if verbose:
        print_status(f"Applied 6σ MAD outlier cleaning: removed {n_outliers}/{len(df)} outliers", "INFO")

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

    if verbose:
        print_status(f"[DATA] Dataset: N={n} observations", "INFO")
        print_status(f"[DATA] Residual mean: {np.mean(y):.6e} m", "INFO")
        print_status(f"[DATA] Residual std:  {np.std(y):.6e} m", "INFO")
        print_status(f"[DATA] Mean uncertainty: {np.mean(y_err):.6e} m", "INFO")
        print_status(f"[DATA] cos(elongation) mean: {np.mean(x):.6f}", "INFO")
        print_status(f"[DATA] cos(elongation) std:  {np.std(x):.6f}", "INFO")

    # MCMC Configuration
    # Compute initial position from fresh OLS on this data — do NOT hardcode
    # the known detection value.  Starting walkers from a data-derived prior is
    # essential so that the Gelman-Rubin R̂ convergence diagnostic is a genuine
    # test of mixing rather than trivially satisfied by pre-positioning walkers.
    X_ols = np.column_stack([x, np.ones_like(x)])
    ols_coeffs, _, _, _ = np.linalg.lstsq(X_ols, y, rcond=None)
    eta_init = ols_coeffs[0] / ETA_SCALE_FACTOR  # η = A / 13
    intercept_init = ols_coeffs[1]
    initial = np.array([eta_init, intercept_init])

    n_walkers = TEP_CONFIG["MCMC_STANDARD_WALKERS"]
    n_steps = TEP_CONFIG["MCMC_STANDARD_STEPS"]
    burn_in = TEP_CONFIG["MCMC_BURN_IN"]
    thin = 20

    if verbose:
        print_status(
            f"  [CALC] OLS-derived initial position: η={eta_init:.6e}, b={intercept_init:.6e}",
            "CALC")

    np.random.seed(42)
    pos = initial + 1e-4 * np.random.randn(n_walkers, 2)

    if verbose:
        print_status("", "INFO")
        print_status("MCMC CONFIGURATION", "PROCESS")
        print_status(f"  [CALC] Number of walkers: {n_walkers}", "CALC")
        print_status(f"  [CALC] Total steps: {n_steps}", "CALC")
        print_status(f"  [CALC] Burn-in steps: {burn_in}", "CALC")
        print_status(f"  [CALC] Thinning factor: {thin}", "CALC")
        print_status(
            f"  [CALC] Initial position: η={initial[0]:.6e}, b={initial[1]:.6e}", "CALC")
        print_status(
            "  [CALC] Prior bounds: η∈[-0.01, 0.01], b∈[-0.1, 0.1]", "CALC")
        print_status(
            f"  [CALC] Expected samples after burn-in/thin: {(n_steps - burn_in) * n_walkers // thin}", "CALC")

    sampler = emcee.EnsembleSampler(
        n_walkers, 2, log_probability, args=(x, y, y_err))

    print_status(
        f"Running MCMC with {n_walkers} walkers for {n_steps} steps...", "INFO")
    sampler.run_mcmc(pos, n_steps, progress=verbose)

    # Get chain for diagnostics
    chain = sampler.get_chain()

    if verbose:
        print_status("", "INFO")
        print_status("MCMC CHAIN DIAGNOSTICS", "PROCESS")

        # Acceptance fraction
        accept_frac = sampler.acceptance_fraction
        print_status("  [CALC] Acceptance fraction per walker:", "CALC")
        for i in range(min(n_walkers, 5)):  # Show first 5
            print_status(
                f"  [CALC]    Walker {i+1:2d}: {accept_frac[i]:.3f}", "CALC")
        print_status(
            f"  [CALC] Mean acceptance: {np.mean(accept_frac):.3f} ± {np.std(accept_frac):.3f}", "CALC")
        print_status(
            "  [CALC] Acceptable range: 0.2-0.5 (ensemble sampler)", "CALC")

        # Gelman-Rubin convergence diagnostic
        R_hat = compute_gelman_rubin(chain)
        print_status(
            "  [CALC] Gelman-Rubin R̂ (convergence diagnostic):", "CALC")
        print_status(f"  [CALC]    R̂_η = {R_hat[0]:.4f} {'✓ converged' if R_hat[0] < 1.1 else '✗ not converged'}",
                     "SUCCESS" if R_hat[0] < 1.1 else "WARNING")
        print_status(f"  [CALC]    R̂_b = {R_hat[1]:.4f} {'✓ converged' if R_hat[1] < 1.1 else '✗ not converged'}",
                     "SUCCESS" if R_hat[1] < 1.1 else "WARNING")
        print_status("  [CALC]    Target: R̂ < 1.1 for convergence", "CALC")

        # Chain statistics by phase
        chain_burned = chain[burn_in:, :, :]
        print_status(
            f"  [CALC] Post-burn-in chain shape: {chain_burned.shape}", "CALC")

        # Evolution of mean by step (trace)
        step_means_eta = np.mean(chain[:, :, 0], axis=1)
        step_stds_eta = np.std(chain[:, :, 0], axis=1)
        print_status("  [CALC] Chain evolution (η):", "CALC")
        print_status(
            f"  [CALC]    Step 0:   mean={step_means_eta[0]:.6e}, std={step_stds_eta[0]:.6e}", "CALC")
        print_status(
            f"  [CALC]    Step 100: mean={step_means_eta[100]:.6e}, std={step_stds_eta[100]:.6e}", "CALC")
        print_status(
            f"  [CALC]    Step 500: mean={step_means_eta[500]:.6e}, std={step_stds_eta[500]:.6e}", "CALC")
        print_status(
            f"  [CALC]    Step 999: mean={step_means_eta[999]:.6e}, std={step_stds_eta[999]:.6e}", "CALC")

    # Get flattened samples after burn-in and thinning
    flat_samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    n_effective = len(flat_samples)

    if verbose:
        print_status("", "INFO")
        print_status("POSTERIOR ANALYSIS", "PROCESS")
        print_status(f"  [CALC] Effective samples: {n_effective}", "CALC")
        print_status(
            f"  [CALC] Thinning reduced samples by factor of {(n_steps - burn_in) * n_walkers / n_effective:.1f}", "CALC")

    # Calculate credible intervals
    eta_samples = flat_samples[:, 0]
    b_samples = flat_samples[:, 1]

    mean_eta = np.mean(eta_samples)
    std_eta = np.std(eta_samples)
    ci_95 = np.percentile(eta_samples, [2.5, 97.5]).tolist()
    ci_68 = np.percentile(eta_samples, [16, 84]).tolist()

    mean_b = np.mean(b_samples)
    std_b = np.std(b_samples)

    if verbose:
        print_status("  [CALC] η posterior statistics:", "CALC")
        print_status(f"  [CALC]    Mean:   {mean_eta:.6e}", "CALC")
        print_status(f"  [CALC]    Std:    {std_eta:.6e}", "CALC")
        print_status(
            f"  [CALC]    Median: {np.median(eta_samples):.6e}", "CALC")
        print_status(
            f"  [CALC]    68% CI: [{ci_68[0]:.6e}, {ci_68[1]:.6e}]", "CALC")
        print_status(
            f"  [CALC]    95% CI: [{ci_95[0]:.6e}, {ci_95[1]:.6e}]", "CALC")
        includes_zero = ci_95[0] <= 0 <= ci_95[1]
        print_status(f"  [CALC]    Includes η=0? {'NO' if not includes_zero else 'YES'}",
                     "SUCCESS" if not includes_zero else "WARNING")

        print_status("  [CALC] b (intercept) posterior statistics:", "CALC")
        print_status(f"  [CALC]    Mean:   {mean_b:.6e} m", "CALC")
        print_status(f"  [CALC]    Std:    {std_b:.6e} m", "CALC")

        # Posterior percentiles for η
        percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
        p_vals = np.percentile(eta_samples, percentiles)
        print_status("  [CALC] η posterior percentiles:", "CALC")
        for p, v in zip(percentiles, p_vals):
            print_status(f"  [CALC]    {p:2d}th: {v:.6e}", "CALC")

    # Bayes Factor Estimate via Savage-Dickey Density Ratio
    # BF = p(η=0|data) / p(η=0) where p(η=0|data) is from KDE at zero
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(eta_samples)

    # Prior density at zero (uniform prior width 0.02 => density = 1/0.02 = 50)
    # Prior range [-0.01, 0.01] for η based on physical constraints:
    # - Nordtvedt effect in GR: η ≈ 10^-13 (negligible)
    # - Alternative theories: η up to 10^-4-10^-3 (e.g., scalar-tensor)
    # - Chosen range [-0.01, 0.01] spans 2 orders of magnitude beyond typical alternative theory predictions,
    #   making it a weakly informative prior that does not constrain the posterior artificially
    prior_density = 1.0 / 0.02  # uniform over [-0.01, 0.01]

    # Posterior density at zero
    posterior_density_at_zero = kde.evaluate(0.0)[0]

    # Savage-Dickey Bayes Factor
    bayes_factor_sd = prior_density / \
        posterior_density_at_zero if posterior_density_at_zero > 0 else np.inf

    # KDE Bandwidth Sensitivity Check
    # Verify that the Savage-Dickey BF is robust to bandwidth choice.
    # Test Scott (default), Silverman, narrow (0.5x), and wide (2.0x).
    bf_sensitivity = {}
    for bw_label, bw_method in [
        ("scott_default", "scott"),
        ("silverman", "silverman"),
        ("narrow_0.5x", 0.5),
        ("wide_2.0x", 2.0),
    ]:
        kde_test = gaussian_kde(eta_samples, bw_method=bw_method)
        post_d_zero = kde_test.evaluate(0.0)[0]
        bf_test = prior_density / post_d_zero if post_d_zero > 0 else np.inf
        bf_sensitivity[bw_label] = float(bf_test)

    bf_sd_values = list(bf_sensitivity.values())
    bf_sd_min = min(bf_sd_values)
    bf_sd_max = max(bf_sd_values)
    bf_sd_range_ratio = bf_sd_max / bf_sd_min if bf_sd_min > 0 else np.inf

    if verbose:
        print_status("  [CALC] Savage-Dickey Bandwidth Sensitivity:", "CALC")
        for label, val in bf_sensitivity.items():
            print_status(f"  [CALC]    {label}: {val:.2e}", "CALC")
        print_status(
            f"  [CALC]    Range: {bf_sd_min:.2e} – {bf_sd_max:.2e} (ratio {bf_sd_range_ratio:.1f}x)",
            "CALC")
        if bf_sd_range_ratio < 10:
            print_status("  [CALC]    Verdict: Robust to bandwidth choice", "SUCCESS")
        else:
            print_status("  [CALC]    Verdict: Sensitive to bandwidth choice", "WARNING")

    # Alternative: BIC approximation
    rss_1 = np.sum((y - (mean_eta * 13.0 * x + mean_b))**2)
    rss_0 = np.sum((y - np.mean(y))**2)
    bic_1 = n * np.log(rss_1/n) + 2 * np.log(n)
    bic_0 = n * np.log(rss_0/n) + 1 * np.log(n)
    delta_bic = bic_0 - bic_1
    bayes_factor_bic = np.exp(0.5 * delta_bic)

    if verbose:
        print_status("", "INFO")
        print_status("BAYES FACTOR CALCULATION", "PROCESS")
        print_status("  [CALC] Savage-Dickey Density Ratio:", "CALC")
        print_status(
            f"  [CALC]    Prior density at η=0:    {prior_density:.2e}", "CALC")
        print_status(
            f"  [CALC]    Posterior KDE at η=0:   {posterior_density_at_zero:.2e}", "CALC")
        print_status(
            f"  [CALC]    BF_SD = prior/post =    {bayes_factor_sd:.2e}", "CALC")

        print_status("  [CALC] BIC Approximation:", "CALC")
        print_status(
            f"  [CALC]    RSS (null model):       {rss_0:.4f}", "CALC")
        print_status(
            f"  [CALC]    RSS (TEP model):        {rss_1:.4f}", "CALC")
        print_status(
            f"  [CALC]    BIC_0 (null):           {bic_0:.2f}", "CALC")
        print_status(
            f"  [CALC]    BIC_1 (TEP):            {bic_1:.2f}", "CALC")
        print_status(
            f"  [CALC]    ΔBIC = BIC_0 - BIC_1:   {delta_bic:.2f}", "CALC")
        print_status(
            f"  [CALC]    BF_BIC = exp(0.5×ΔBIC): {bayes_factor_bic:.2e}", "CALC")

        # Jeffreys scale interpretation
        log_bf = np.log10(bayes_factor_bic)
        if log_bf > 2:
            evidence = "Decisive"
        elif log_bf > 1.5:
            evidence = "Very Strong"
        elif log_bf > 1:
            evidence = "Strong"
        elif log_bf > 0.5:
            evidence = "Substantial"
        else:
            evidence = "Weak/None"
    
    # Bayes factor interpretation thresholds from Kass & Raftery (1995):
    # log10(BF) > 1: Strong evidence, > 0.5: Substantial evidence, < 0.5: Weak/None

        print_status(f"  [CALC] log₁₀(BF) = {log_bf:.2f}", "CALC")
        print_status(
            f"  [CALC] Evidence strength: {evidence}", "SUCCESS" if log_bf > 1 else "INFO")
        print_status(
            f"  [CALC] P(TEP|data) with equal priors: {bayes_factor_bic/(1+bayes_factor_bic):.6f}", "CALC")

    results = {
        "posterior_mean_eta": float(mean_eta),
        "posterior_std_eta": float(std_eta),
        "posterior_mean_b": float(mean_b),
        "posterior_std_b": float(std_b),
        "credible_interval_95": ci_95,
        "credible_interval_68": ci_68,
        "bayes_factor_savage_dickey": float(bayes_factor_sd),
        "bayes_factor_bic": float(bayes_factor_bic),
        "bayes_factor_sensitivity": {
            "bandwidth_methods": bf_sensitivity,
            "min": float(bf_sd_min),
            "max": float(bf_sd_max),
            "range_ratio": float(bf_sd_range_ratio),
            "robust": bool(bf_sd_range_ratio < 10) if np.isfinite(bf_sd_range_ratio) else None
        },
        "n_samples": n_effective,
        "gelman_rubin_eta": float(R_hat[0]) if 'R_hat' in locals() else None,
        "gelman_rubin_b": float(R_hat[1]) if 'R_hat' in locals() else None,
        "mean_acceptance_fraction": float(np.mean(accept_frac)) if 'accept_frac' in locals() else None,
        "outlier_cleaning": {
            "method": "6σ MAD (detect_outliers_sigma)",
            "sigma_threshold": 6.0,
            "n_outliers_removed": n_outliers,
            "n_original": len(df),
            "n_cleaned": len(df_clean)
        }
    }

    if verbose:
        print_status("", "INFO")
        print_status("="*60, "INFO")
        print_status("BAYESIAN ANALYSIS SUMMARY", "TITLE")
        print_status("="*60, "INFO")
        print_status(
            f"  Posterior η: {results['posterior_mean_eta']:.6e} ± {results['posterior_std_eta']:.6e}", "CALC")
        print_status(
            f"  95% CI:      [{results['credible_interval_95'][0]:.6e}, {results['credible_interval_95'][1]:.6e}]", "CALC")
        print_status(
            f"  Bayes Factor: {results['bayes_factor_bic']:.2e} ({evidence} evidence)", "SUCCESS" if log_bf > 1 else "CALC")

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 016: Bayesian MCMC Analysis")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_016", str(
        log_dir / "step_016_bayesian_analysis.log"))
    set_step_logger(logger)

    print_status("Starting Bayesian MCMC Analysis...", "TITLE")

    # run_bayesian_analysis loads data internally
    summary = run_bayesian_analysis(verbose=True)

    results = {
        "step_id": "step_016",
        "bayesian_summary": summary,
        "status": "PASS" if summary["bayes_factor_bic"] > 10 else "WARNING"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_016_bayesian_analysis")
    print_status("Bayesian Analysis Complete.", "SUCCESS")

