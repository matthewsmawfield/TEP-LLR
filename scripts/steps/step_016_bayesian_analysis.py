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
from scripts.utils.numerics import stable_lstsq
import pandas as pd
import emcee
from scipy import stats
from scipy.stats import gaussian_kde

from scripts.utils.bayesian_evidence import (
    ETA_PRIOR_SPECS,
    build_prior_sensitivity_table,
)
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.statistical_utils import detect_outliers_sigma

TEP_CONFIG = get_config()


def gelman_rubin(chains):
    """Compute Gelman-Rubin statistic for convergence diagnostics."""
    # emcee stores (n_steps, n_walkers, n_params); statistics need (n_walkers, n_steps, n_params)
    if chains.shape[0] > chains.shape[1]:
        chains = np.transpose(chains, (1, 0, 2))
    n_walkers, n_steps, _n_params = chains.shape

    W = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    chain_means = np.mean(chains, axis=1)
    B = n_steps * np.var(chain_means, axis=0, ddof=1)

    # Estimated variance
    V_hat = ((n_steps - 1) / n_steps) * W + (1 / n_steps) * B

    # Gelman-Rubin statistic
    R_hat = np.sqrt(V_hat / W)

    return R_hat


# Alias for compatibility
compute_gelman_rubin = gelman_rubin


def log_probability(theta, x, y, y_err, eta_lo=-0.01, eta_hi=0.01):
    """Log posterior for synodic-only model with uniform box prior on (η, b)."""
    eta, intercept = theta
    if not (eta_lo <= eta <= eta_hi and -0.1 <= intercept <= 0.1):
        return -np.inf
    model = eta * ETA_SCALE_FACTOR * x + intercept
    residuals = y - model
    chi2 = np.sum((residuals / y_err) ** 2)
    return -0.5 * chi2


def _run_mcmc_chain(
    x,
    y,
    y_err,
    eta_lo,
    eta_hi,
    n_walkers,
    n_steps,
    burn_in,
    thin,
    verbose,
    label,
):
    """Run emcee for one documented η prior; return flat η samples and diagnostics."""
    X_ols = np.column_stack([x, np.ones_like(x)])
    ols_coeffs, _, _, _ = stable_lstsq(X_ols, y)
    eta_init = ols_coeffs[0] / ETA_SCALE_FACTOR
    intercept_init = ols_coeffs[1]
    initial = np.array([eta_init, intercept_init])

    np.random.seed(TEP_CONFIG.get("RANDOM_SEED", 42))
    pos = initial + 1e-4 * np.random.randn(n_walkers, 2)

    sampler = emcee.EnsembleSampler(
        n_walkers,
        2,
        log_probability,
        args=(x, y, y_err, eta_lo, eta_hi),
    )
    if verbose:
        print_status(
            f"MCMC [{label}]: η∈[{eta_lo}, {eta_hi}], {n_walkers} walkers × {n_steps} steps",
            "INFO",
        )
    sampler.run_mcmc(pos, n_steps, progress=verbose)

    chain = sampler.get_chain()
    R_hat = compute_gelman_rubin(chain)
    flat_samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    return {
        "flat_samples": flat_samples,
        "eta_samples": flat_samples[:, 0],
        "R_hat": R_hat,
        "acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
    }

def run_bayesian_analysis(verbose=False):
    if verbose:
        print_status("="*60, "INFO")
        print_status("BAYESIAN MCMC ANALYSIS - DETAILED TRACE", "TITLE")
        print_status("="*60, "INFO")

    input_path = PROJECT_ROOT / "data" / "processed" / \
        "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)

    # Apply 6σ-equivalent (MAD-based) outlier cleaning for consistency with other steps (CRITICAL FIX for standardization)
    outlier_mask = detect_outliers_sigma(df['residual_m'].values, sigma_threshold=6.0)
    n_outliers = int(np.sum(outlier_mask))
    df_clean = df[~outlier_mask]  # PERFORMANCE FIX: Removed unnecessary .copy()
    if verbose:
        print_status(f"Applied 6σ-equivalent (MAD-based) outlier cleaning: removed {n_outliers}/{len(df)} outliers", "INFO")

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

    n_walkers = TEP_CONFIG["MCMC_STANDARD_WALKERS"]
    n_steps_ref = TEP_CONFIG["MCMC_STANDARD_STEPS"]
    n_steps_sensitivity = max(2000, n_steps_ref // 2)
    burn_in = TEP_CONFIG["MCMC_BURN_IN"]
    thin = 20

    print_status(
        f"Running synodic MCMC: reference {n_steps_ref} steps; "
        f"sensitivity priors {n_steps_sensitivity} steps each",
        "INFO",
    )

    posterior_by_prior = {}
    flat_by_prior = {}
    mcmc_runs = {}
    for spec in ETA_PRIOR_SPECS:
        steps = n_steps_ref if spec.get("reference") else n_steps_sensitivity
        run = _run_mcmc_chain(
            x,
            y,
            y_err,
            spec["eta_lo"],
            spec["eta_hi"],
            n_walkers,
            steps,
            burn_in,
            thin,
            verbose,
            spec["prior_id"],
        )
        posterior_by_prior[spec["prior_id"]] = run["eta_samples"]
        flat_by_prior[spec["prior_id"]] = run["flat_samples"]
        mcmc_runs[spec["prior_id"]] = {
            "n_steps": steps,
            "gelman_rubin_eta": float(run["R_hat"][0]),
            "gelman_rubin_b": float(run["R_hat"][1]),
            "mean_acceptance_fraction": run["acceptance_fraction"],
            "n_posterior_samples": int(len(run["eta_samples"])),
        }

    ref_id = "uniform_reference"
    flat_samples = flat_by_prior[ref_id]
    n_effective = len(flat_samples)
    R_hat = np.array(
        [
            mcmc_runs[ref_id]["gelman_rubin_eta"],
            mcmc_runs[ref_id]["gelman_rubin_b"],
        ]
    )
    accept_frac_mean = mcmc_runs[ref_id]["mean_acceptance_fraction"]

    if verbose:
        print_status("", "INFO")
        print_status("POSTERIOR ANALYSIS", "PROCESS")
        print_status(f"  [CALC] Effective samples: {n_effective}", "CALC")
        print_status(
            f"  [CALC] Thinning reduced samples by factor of {(n_steps_ref - burn_in) * n_walkers / n_effective:.1f}", "CALC")

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

    prior_sensitivity_table = build_prior_sensitivity_table(
        posterior_by_prior, ETA_PRIOR_SPECS
    )
    ref_row = next(r for r in prior_sensitivity_table if r["prior_id"] == ref_id)
    bf_sensitivity = ref_row["bandwidth_bf"]
    prior_density = ref_row["prior_density_at_eta_zero"]
    bayes_factor_sd = float(ref_row["bandwidth_bf"]["scott_default"])
    bf_sd_values = list(bf_sensitivity.values())
    bf_sd_min = ref_row["bf_min"]
    bf_sd_max = ref_row["bf_max"]
    bf_sd_range_ratio = ref_row["bf_range_ratio"]
    bf_sd_gmean = ref_row["bf_geometric_mean"]
    posterior_density_at_zero = float(
        gaussian_kde(eta_samples).evaluate(0.0)[0]
    )

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

    # Alternative: BIC approximation (robust, no bandwidth ambiguity)
    rss_1 = np.sum((y - (mean_eta * ETA_SCALE_FACTOR * x + mean_b))**2)
    rss_0 = np.sum((y - np.mean(y))**2)
    bic_1 = n * np.log(rss_1/n) + 2 * np.log(n)
    bic_0 = n * np.log(rss_0/n) + 1 * np.log(n)
    delta_bic = bic_0 - bic_1
    bayes_factor_bic = np.exp(0.5 * delta_bic)

    # Primary Bayes Factor: BIC approximation (robust, no bandwidth ambiguity)
    # Savage-Dickey is reported as sensitivity analysis only.
    bayes_factor_primary = float(bayes_factor_bic)

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
        "bayes_factor_primary": bayes_factor_primary,
        "bayes_factor_primary_method": "BIC approximation (robust to bandwidth ambiguity)",
        "bayes_factor_bic": float(bayes_factor_bic),
        "bayes_factor_savage_dickey": float(bayes_factor_sd),
        "bayes_factor_savage_dickey_geometric_mean": float(bf_sd_gmean),
        "bayes_factor_sensitivity": {
            "note": "Savage-Dickey is bandwidth- and prior-sensitive; BIC is the primary robust estimator.",
            "bandwidth_methods": bf_sensitivity,
            "min": float(bf_sd_min),
            "max": float(bf_sd_max),
            "geometric_mean": float(bf_sd_gmean),
            "range_ratio": float(bf_sd_range_ratio),
            "robust": bool(bf_sd_range_ratio < 10) if np.isfinite(bf_sd_range_ratio) else None,
        },
        "prior_sensitivity_table": prior_sensitivity_table,
        "prior_specifications": ETA_PRIOR_SPECS,
        "mcmc_runs_by_prior": mcmc_runs,
        "n_samples": n_effective,
        "gelman_rubin_eta": float(R_hat[0]),
        "gelman_rubin_b": float(R_hat[1]),
        "mean_acceptance_fraction": float(accept_frac_mean),
        "outlier_cleaning": {
            "method": "6σ-equivalent (MAD-based; detect_outliers_sigma)",
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

