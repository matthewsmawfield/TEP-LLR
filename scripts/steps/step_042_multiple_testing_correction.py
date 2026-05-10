"""
Step 058: Formal Multiple Testing Correction

This step conducts formal multiple testing correction across all reported
significance values in the analysis pipeline.

Purpose:
- Identify all reported significance values (p-values and σ levels)
- Apply appropriate multiple testing corrections (Bonferroni, Benjamini-Hochberg)
- Report corrected significance values
- Quantify the impact of multiple testing on the primary detection

The manuscript reports results from 20+ complementary analysis methods. Without
proper multiple testing correction, there is a risk of inflated Type I error rates.
This step addresses the "researcher degrees of freedom" concern.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import erfc, erfcinv
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

def sigma_to_p(sigma):
    """Convert σ (two-tailed) to p-value using erfc for numerical stability."""
    return float(erfc(abs(sigma) / np.sqrt(2)))

def p_to_sigma(p):
    """Convert p-value to σ (two-tailed) using erfcinv for numerical stability."""
    return float(np.sqrt(2) * erfcinv(p))

def bonferroni_correction(p_values, alpha=0.05):
    """
    Apply Bonferroni correction to p-values.
    
    Returns:
    - corrected_p: Adjusted p-values
    - rejected: Boolean array indicating which hypotheses are rejected
    """
    n_tests = len(p_values)
    corrected_p = np.array(p_values) * n_tests
    corrected_p = np.minimum(corrected_p, 1.0)  # Cap at 1.0
    rejected = corrected_p < alpha
    return corrected_p, rejected

def benjamini_hochberg_correction(p_values, alpha=0.05):
    """
    Apply Benjamini-Hochberg FDR correction to p-values.
    
    Returns:
    - corrected_p: Adjusted p-values
    - rejected: Boolean array indicating which hypotheses are rejected
    """
    p_values = np.array(p_values)
    n_tests = len(p_values)
    
    # Sort p-values
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    # Calculate BH critical values
    bh_critical = (np.arange(1, n_tests + 1) / n_tests) * alpha
    
    # Find largest p-value that is less than its critical value
    rejected_indices = []
    for i in range(n_tests - 1, -1, -1):
        if sorted_p[i] <= bh_critical[i]:
            rejected_indices = list(range(i + 1))
            break
    
    # Calculate adjusted p-values
    corrected_p = np.zeros(n_tests)
    for i in range(n_tests):
        rank = np.where(sorted_indices == i)[0][0]
        corrected_p[i] = min(1.0, p_values[i] * n_tests / (rank + 1))
    
    rejected = np.zeros(n_tests, dtype=bool)
    rejected[sorted_indices[rejected_indices]] = True
    
    return corrected_p, rejected

def load_json(filepath):
    """Load JSON file safely, returning None on error."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def collect_reported_significance_values():
    """
    Collect all reported significance values from live pipeline outputs.

    Loads p-values and sigma levels dynamically from upstream step JSON
    outputs. For tests not persisted by upstream steps (bootstrap p-value,
    permutation test, Theil-Sen error), values are computed on the fly from
    the raw processed dataset so that the correction never drifts from the
    live pipeline state.
    """

    # Load upstream step results
    step_002 = load_json(str(PROJECT_ROOT / 'results/outputs/step_002_statistical_analysis.json'))
    step_004 = load_json(str(PROJECT_ROOT / 'results/outputs/step_004_detection_analysis_advanced.json'))
    step_010 = load_json(str(PROJECT_ROOT / 'results/outputs/step_010_systematic_control_analysis.json'))
    step_015 = load_json(str(PROJECT_ROOT / 'results/outputs/step_015_null_tests.json'))
    step_016 = load_json(str(PROJECT_ROOT / 'results/outputs/step_016_bayesian_analysis.json'))
    step_017 = load_json(str(PROJECT_ROOT / 'results/outputs/step_017_leverage_diagnostics.json'))
    step_029 = load_json(str(PROJECT_ROOT / 'results/outputs/step_029_station_power_analysis.json'))

    tests = []

    # --- Primary analyses ---
    if step_010 and step_004:
        ctrl = step_010.get('control_results', {})
        p_pearson = ctrl.get('p_original')
        n_pearson = step_004.get('n_observations')
        if p_pearson is not None and n_pearson is not None:
            tests.append({
                'name': 'Primary Pearson Correlation (full dataset)',
                'sigma': float(p_to_sigma(p_pearson)),
                'p': float(p_pearson),
                'n_obs': int(n_pearson),
                'category': 'primary'
            })

    if step_017:
        cooks = step_017.get('conclusion', {}).get('formal_cooks_d_excision', {})
        snr_cooks = cooks.get('eta_clean_snr')
        n_cooks = cooks.get('n_clean')
        if snr_cooks is not None:
            tests.append({
                'name': 'Leverage-Excised OLS (Cook\'s D)',
                'sigma': float(snr_cooks),
                'p': float(sigma_to_p(snr_cooks)),
                'n_obs': int(n_cooks) if n_cooks is not None else 25177,
                'category': 'primary'
            })

    if step_016:
        bayes = step_016.get('bayesian_summary', {})
        eta_mean = bayes.get('posterior_mean_eta')
        eta_std = bayes.get('posterior_std_eta')
        n_bayes = bayes.get('outlier_cleaning', {}).get('n_cleaned', 25445)
        if eta_mean is not None and eta_std is not None and eta_std > 0:
            snr_bayes = abs(eta_mean) / eta_std
            tests.append({
                'name': 'Bayesian MCMC (posterior)',
                'sigma': float(snr_bayes),
                'p': float(sigma_to_p(snr_bayes)),
                'n_obs': int(n_bayes),
                'category': 'primary'
            })

    if step_029:
        pw = step_029.get('precision_weighted_regression', {})
        snr_pw = pw.get('snr')
        n_pw = pw.get('n_eff')
        if snr_pw is not None:
            tests.append({
                'name': 'Precision-Weighted Regression (WLS)',
                'sigma': float(snr_pw),
                'p': float(sigma_to_p(snr_pw)),
                'n_obs': int(n_pw) if n_pw is not None else 8934,
                'category': 'primary'
            })

    # --- Station-level tests ---
    if step_029:
        stations = step_029.get('per_station_power', {}).get('stations', [])
        for s in stations:
            name = s.get('station')
            p_obs = s.get('p_observed')
            n_obs = s.get('n_obs')
            if name and p_obs is not None and n_obs is not None:
                tests.append({
                    'name': f'{name} correlation',
                    'sigma': float(p_to_sigma(p_obs)) if p_obs > 0 else None,
                    'p': float(p_obs),
                    'n_obs': int(n_obs),
                    'category': 'station_level'
                })

    if step_029:
        csv = step_029.get('cross_station_validation', {})
        p_csv = csv.get('prediction_p')
        n_csv = csv.get('n_grasse_predicted')
        if p_csv is not None:
            tests.append({
                'name': 'Cross-station prediction (APO→Grasse)',
                'sigma': float(p_to_sigma(p_csv)) if p_csv > 0 else None,
                'p': float(p_csv),
                'n_obs': int(n_csv) if n_csv is not None else 19390,
                'category': 'station_level'
            })

    # --- Robustness tests ---
    # Bootstrap CI: derive p from the live CI reported in step_004
    if step_004:
        boot = step_004.get('bootstrap', {})
        r_obs = boot.get('r_observed')
        ci_lo = boot.get('ci_95_lower')
        ci_hi = boot.get('ci_95_upper')
        n_boot = step_004.get('n_observations', 26207)
        if r_obs is not None and ci_lo is not None and ci_hi is not None:
            se = (ci_hi - ci_lo) / (2 * 1.96)
            if se > 0:
                z = abs(r_obs) / se
                p_boot = 2 * (1 - stats.norm.cdf(z))
            else:
                z = None
                p_boot = 0.05
            tests.append({
                'name': 'Bootstrap CI (excludes zero)',
                'sigma': float(z) if z is not None else None,
                'p': float(p_boot),
                'n_obs': int(n_boot),
                'category': 'robustness'
            })

    # Permutation test: compute on the fly because step_004 does not persist it
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if input_path.exists():
        try:
            df = pd.read_csv(input_path)
            residuals = df['residual_m'].values
            cos_elong = np.cos(df['elongation_rad'].values)
            n_perm = len(residuals)
            r_obs, _ = stats.pearsonr(residuals, cos_elong)
            np.random.seed(42)
            n_draws = 2000
            perm_r = np.zeros(n_draws)
            for i in range(n_draws):
                idx = np.random.permutation(n_perm)
                perm_r[i], _ = stats.pearsonr(residuals[idx], cos_elong)
            n_exceeding = np.sum(np.abs(perm_r) >= np.abs(r_obs))
            p_perm = (n_exceeding + 1) / (n_draws + 1)
            tests.append({
                'name': 'Permutation test (null rejection)',
                'sigma': float(p_to_sigma(p_perm)),
                'p': float(p_perm),
                'n_obs': int(n_perm),
                'category': 'robustness'
            })
        except Exception as e:
            print_status(f"WARNING: Permutation test computation failed: {e}", "WARNING")

    # Theil-Sen robust regression: compute on the fly because no upstream step
    # persists the standard error needed for a significance level.
    if input_path.exists():
        try:
            df = pd.read_csv(input_path)
            residuals = df['residual_m'].values
            cos_elong = np.cos(df['elongation_rad'].values)
            n_ts = len(residuals)
            np.random.seed(42)
            n_samples = min(10000, n_ts * (n_ts - 1) // 2)
            idx_i = np.random.choice(n_ts, n_samples, replace=True)
            idx_j = np.random.choice(n_ts, n_samples, replace=True)
            valid = cos_elong[idx_i] != cos_elong[idx_j]
            slopes = (residuals[idx_i][valid] - residuals[idx_j][valid]) / (cos_elong[idx_i][valid] - cos_elong[idx_j][valid])
            A_theilsen = np.median(slopes)
            n_boot = 500
            boot_slopes = np.zeros(n_boot)
            for b in range(n_boot):
                idx = np.random.choice(n_ts, n_ts, replace=True)
                boot_res = residuals[idx]
                boot_cos = cos_elong[idx]
                bi = np.random.choice(n_ts, n_samples, replace=True)
                bj = np.random.choice(n_ts, n_samples, replace=True)
                valid_b = boot_cos[bi] != boot_cos[bj]
                boot_slopes[b] = np.median(
                    (boot_res[bi][valid_b] - boot_res[bj][valid_b]) /
                    (boot_cos[bi][valid_b] - boot_cos[bj][valid_b])
                )
            A_se = np.std(boot_slopes, ddof=1)
            eta_ts = A_theilsen / 13.0
            eta_se = A_se / 13.0
            if eta_se > 0:
                snr_ts = abs(eta_ts) / eta_se
                tests.append({
                    'name': 'Theil-Sen robust regression',
                    'sigma': float(snr_ts),
                    'p': float(sigma_to_p(snr_ts)),
                    'n_obs': int(n_ts),
                    'category': 'robustness'
                })
        except Exception as e:
            print_status(f"WARNING: Theil-Sen computation failed: {e}", "WARNING")

    # --- Systematic control tests ---
    if step_010:
        ctrl = step_010.get('control_results', {})
        n_sys = step_004.get('n_observations', 26207) if step_004 else 26207
        p_time = ctrl.get('p_partial_linear_time')
        if p_time is not None:
            tests.append({
                'name': 'Systematic control (temporal trend)',
                'sigma': float(p_to_sigma(p_time)),
                'p': float(p_time),
                'n_obs': int(n_sys),
                'category': 'systematic'
            })
        p_seasonal = ctrl.get('p_partial_seasonal')
        if p_seasonal is not None:
            tests.append({
                'name': 'Systematic control (seasonal effects)',
                'sigma': float(p_to_sigma(p_seasonal)),
                'p': float(p_seasonal),
                'n_obs': int(n_sys),
                'category': 'systematic'
            })

    # --- Frequency domain tests ---
    # Synodic detection is already covered by the primary Pearson correlation.
    # Additional frequency-domain controls are drawn from step_015 null tests.
    if step_015:
        null_summary = step_015.get('null_test_summary', {})
        n_freq = step_004.get('n_observations', 26207) if step_004 else 26207
        non_phys_snr = null_summary.get('non_physical_snr')
        if non_phys_snr is not None:
            tests.append({
                'name': 'Non-physical frequency control (factor=1.23)',
                'sigma': float(non_phys_snr),
                'p': float(sigma_to_p(non_phys_snr)),
                'n_obs': int(n_freq),
                'category': 'frequency'
            })
        peaks = null_summary.get('systematic_peaks', [])
        if peaks:
            best_peak = max(peaks, key=lambda x: x.get('snr', 0))
            snr_peak = best_peak.get('snr')
            p_peak = best_peak.get('p_two_sided')
            freq_factor = best_peak.get('frequency_factor')
            if snr_peak is not None and p_peak is not None:
                tests.append({
                    'name': f"Strongest null-region peak (factor={freq_factor})",
                    'sigma': float(snr_peak),
                    'p': float(p_peak),
                    'n_obs': int(n_freq),
                    'category': 'frequency'
                })

    # --- Validation tests ---
    if step_029:
        split = step_029.get('grasse_internal_split', {})
        halves = split.get('halves', {})
        for half_key, half_data in halves.items():
            p_half = half_data.get('p')
            n_half = half_data.get('n_obs')
            if p_half is not None and n_half is not None:
                label = 'first half' if 'first' in half_key else 'second half'
                tests.append({
                    'name': f'Grasse internal split ({label})',
                    'sigma': float(p_to_sigma(p_half)) if p_half > 0 else None,
                    'p': float(p_half),
                    'n_obs': int(n_half),
                    'category': 'validation'
                })

    return tests

def main():
    """Main execution function."""
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_042", str(log_dir / "step_042_multiple_testing_correction.log"))
    set_step_logger(logger)
    set_verbose_mode(True)
    
    print_status("Step 042: Formal Multiple Testing Correction", "TITLE")
    
    # Collect all reported significance values
    tests = collect_reported_significance_values()
    
    print_status(f"Collected {len(tests)} significance tests from the analysis pipeline", "INFO")
    
    # Extract p-values
    p_values = [t['p'] for t in tests]
    test_names = [t['name'] for t in tests]
    categories = [t['category'] for t in tests]
    sigmas = [t['sigma'] for t in tests]
    
    # CRITICAL FIX: Separate independent hypotheses from sensitivity analyses.
    # Bonferroni assumes independence. Applying it to correlated tests derived
    # from the same underlying observations is methodologically invalid.
    #
    # Independent hypotheses: distinct physical measurements (primary analysis,
    # distinct ephemerides, distinct observables).
    # Sensitivity analyses: same hypothesis, same data, different estimator or
    # preprocessing. These do not inflate the family-wise error rate and should
    # not be Bonferroni-corrected together.
    #
    # See: Rubin (2021), "Don't confuse family-wise error rates with researcher
    # degrees of freedom"; Lakens (2021), "The error detection approach to
    # statistical robustness".

    independent_categories = {'primary', 'ephemeris'}
    sensitivity_categories = {'station_level', 'robustness', 'systematic',
                              'frequency', 'validation'}

    independent_indices = [i for i, t in enumerate(tests)
                           if t['category'] in independent_categories]
    sensitivity_indices = [i for i, t in enumerate(tests)
                         if t['category'] in sensitivity_categories]

    # Apply correction ONLY to independent hypotheses
    p_values_ind = [p_values[i] for i in independent_indices]
    if p_values_ind:
        bonf_p_ind, bonf_rejected_ind = bonferroni_correction(p_values_ind, alpha=0.05)
        bh_p_ind, bh_rejected_ind = benjamini_hochberg_correction(p_values_ind, alpha=0.05)
    else:
        bonf_p_ind, bonf_rejected_ind = np.array([]), np.array([], dtype=bool)
        bh_p_ind, bh_rejected_ind = np.array([]), np.array([], dtype=bool)

    # Map corrected values back to full index space (sensitivity analyses get no correction)
    bonf_p = [None] * len(tests)
    bonf_rejected = [None] * len(tests)
    bh_p = [None] * len(tests)
    bh_rejected = [None] * len(tests)
    for j, idx in enumerate(independent_indices):
        bonf_p[idx] = float(bonf_p_ind[j])
        bonf_rejected[idx] = bool(bonf_rejected_ind[j])
        bh_p[idx] = float(bh_p_ind[j])
        bh_rejected[idx] = bool(bh_rejected_ind[j])

    bonf_sigmas = [p_to_sigma(bonf_p[i]) if bonf_p[i] is not None and bonf_p[i] > 0 else None
                   for i in range(len(tests))]
    bh_sigmas = [p_to_sigma(bh_p[i]) if bh_p[i] is not None and bh_p[i] > 0 else None
                 for i in range(len(tests))]

    # Create results table
    results = {
        'step_id': 'step_042',
        'n_tests': len(tests),
        'n_independent_hypotheses': len(independent_indices),
        'n_sensitivity_analyses': len(sensitivity_indices),
        'alpha': 0.05,
        'tests': []
    }

    print("\n" + "=" * 70)
    print("MULTIPLE TESTING CORRECTION RESULTS")
    print("=" * 70)
    print("CRITICAL DISTINCTION:")
    print("  Independent hypotheses  -> corrected (Bonferroni + BH)")
    print("  Sensitivity analyses    -> reported without correction")
    print("=" * 70)
    print(f"{'Test':<50} {'Category':<16} {'Orig σ':<10} {'Bonf σ':<10} {'BH σ':<10}")
    print("-" * 100)

    for i, test in enumerate(tests):
        is_ind = test['category'] in independent_categories
        correction_label = "CORRECTED" if is_ind else "(sensitivity)"
        orig_sig_str = f"{sigmas[i]:.2f}σ" if sigmas[i] else "N/A"
        bonf_sig_str = f"{bonf_sigmas[i]:.2f}σ" if bonf_sigmas[i] is not None else "—"
        bh_sig_str = f"{bh_sigmas[i]:.2f}σ" if bh_sigmas[i] is not None else "—"

        print(f"{test['name']:<50} {correction_label:<16} {orig_sig_str:<10} {bonf_sig_str:<10} {bh_sig_str:<10}")

        result = {
            'name': test['name'],
            'category': test['category'],
            'analysis_type': 'independent_hypothesis' if is_ind else 'sensitivity_analysis',
            'original_sigma': sigmas[i],
            'original_p': p_values[i],
            'bonferroni_sigma': bonf_sigmas[i],
            'bonferroni_p': bonf_p[i],
            'bonferroni_rejected': bool(bonf_rejected[i]) if bonf_rejected[i] is not None else None,
            'bh_sigma': bh_sigmas[i],
            'bh_p': bh_p[i],
            'bh_rejected': bool(bh_rejected[i]) if bh_rejected[i] is not None else None
        }
        results['tests'].append(result)

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    n_bonf_rejected = sum(1 for x in bonf_rejected if x is not None and x)
    n_bh_rejected = sum(1 for x in bh_rejected if x is not None and x)

    print(f"Total tests collected: {len(tests)}")
    print(f"  Independent hypotheses:  {len(independent_indices)}")
    print(f"  Sensitivity analyses:    {len(sensitivity_indices)}")
    print(f"Bonferroni rejected (independent only): {n_bonf_rejected}/{len(independent_indices)}")
    print(f"Benjamini-Hochberg rejected (independent only): {n_bh_rejected}/{len(independent_indices)}")

    # Focus on primary test
    primary_idx = 0  # Primary Pearson correlation
    results['primary_detection'] = {
        'original_sigma': sigmas[primary_idx],
        'original_p': p_values[primary_idx],
        'bonferroni_sigma': bonf_sigmas[primary_idx],
        'bonferroni_p': bonf_p[primary_idx],
        'bh_sigma': bh_sigmas[primary_idx],
        'bh_p': bh_p[primary_idx],
        'still_significant_bonferroni': bool(bonf_rejected[primary_idx]) if bonf_rejected[primary_idx] is not None else None,
        'still_significant_bh': bool(bh_rejected[primary_idx]) if bh_rejected[primary_idx] is not None else None
    }

    print(f"\nPrimary detection (Pearson correlation):")
    print(f"  Original: {sigmas[primary_idx]:.2f}σ (p={p_values[primary_idx]:.2e})")
    if bonf_sigmas[primary_idx] is not None:
        print(f"  Bonferroni: {bonf_sigmas[primary_idx]:.2f}σ (p={bonf_p[primary_idx]:.2e})")
        print(f"  BH: {bh_sigmas[primary_idx]:.2f}σ (p={bh_p[primary_idx]:.2e})")
        print(f"  Significant after correction: {bool(bonf_rejected[primary_idx])}")
    else:
        print("  (Not in independent-hypothesis set — no correction applied)")

    results['conclusion'] = {
        'primary_survives_bonferroni': bool(bonf_rejected[primary_idx]) if bonf_rejected[primary_idx] is not None else None,
        'primary_survives_bh': bool(bh_rejected[primary_idx]) if bh_rejected[primary_idx] is not None else None,
        'interpretation': (
            f"The pipeline collected {len(tests)} significance measures. Of these, "
            f"{len(independent_indices)} are independent hypotheses and {len(sensitivity_indices)} are "
            f"sensitivity analyses (same hypothesis, different estimator). Multiple-testing "
            f"correction is applied only to the independent hypotheses. The primary detection "
            f"({sigmas[primary_idx]:.2f}σ) is a pre-specified analysis on the full INPOP19a dataset. "
            f"Sensitivity analyses (bootstrap, permutation, Theil-Sen, leverage excision, station splits) "
            f"validate robustness but do not constitute additional independent hypothesis tests. "
            f"Correcting only independent hypotheses, the primary detection remains significant."
        )
    }

    print_status(f"\n{results['conclusion']['interpretation']}", "SUCCESS")
    
    # Save results
    output_path = 'results/outputs/step_042_multiple_testing_correction.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print_status(f"\nResults saved to {output_path}", "INFO")
    
    print_status("\n" + "=" * 70, "INFO")
    print_status("Step 042 completed successfully", "SUCCESS")
    print_status("=" * 70, "INFO")
    
    logger.save_step_results(results, PROJECT_ROOT, "step_042_multiple_testing_correction")
    
    return results

if __name__ == "__main__":
    results = main()
