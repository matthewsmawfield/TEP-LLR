#!/usr/bin/env python3
"""
Statistical utilities for TEP-LLR analysis.

Contains reusable functions for outlier detection, power analysis,
linear regression, and other statistical methods used across the pipeline.
"""

import numpy as np
from scipy import linalg, stats
from typing import Dict

# Import constants from llr_constants module
from scripts.utils.llr_constants import ETA_SCALE_FACTOR, Z_ALPHA_2
from scripts.utils.config import get_config

TEP_CONFIG = get_config()

def robust_regression(y: np.ndarray, X: np.ndarray, weights: np.ndarray = None, 
                      scale_errors_by_birge: bool = True) -> Dict:
    """
    Perform numerically stable weighted linear regression using QR decomposition.
    
    Formula:
    --------
    Solve for β: y = Xβ + ε
    Weighted version: W^{1/2}y = W^{1/2}Xβ + ε'
    Using QR: Rβ = Qᵀ(W^{1/2}y)
    
    Covariance:
    -----------
    Σ_β = (XᵀWX)⁻¹ σ² = R⁻¹(R⁻ᵀ) σ²
    where σ² is the unbiased variance estimate (Birge-scaled if enabled).
    
    Parameters:
    -----------
    y : np.ndarray, shape (n,)
        Dependent variable (residuals) [m]
    X : np.ndarray, shape (n, k)
        Design matrix (including intercept column if needed)
    weights : np.ndarray, shape (n,), optional
        Observation weights (1/σ²). If None, uniform weighting is used.
    scale_errors_by_birge : bool, default True
        If True, scales formal errors by the Birge Ratio max(1.0, sqrt(chi2_red)).
        
    Returns:
    --------
    Dict containing:
        - coefficients: np.ndarray (β)
        - errors: np.ndarray (scaled standard errors)
        - chi2_red: Reduced chi-squared
        - birge_ratio: Birge Ratio (sqrt of chi2_red)
        - condition_number: Matrix condition number κ(R)
        - n_obs: Number of observations
        - dof: Degrees of freedom
        - rss: Residual sum of squares
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, k = X.shape
    dof = n - k
    
    if n <= k:
        return {'coefficients': np.full(k, np.nan), 'errors': np.full(k, np.nan), 
                'chi2_red': np.nan, 'birge_ratio': np.nan, 'condition_number': np.nan}
        
    if weights is None:
        weights = np.ones(n)
    
    # Weighted transformation
    sqrt_w = np.sqrt(weights)
    yw = y * sqrt_w
    Xw = X * sqrt_w[:, np.newaxis]
    
    try:
        # QR decomposition and covariance: suppress benign LAPACK overflow
        # warnings from near-rank-deficient weighted designs (SciPy / NumPy 2.x).
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            Q, R = linalg.qr(Xw, mode='economic')

            # Check condition number
            # Stability threshold: kappa < 1e12 for float64
            s = linalg.svdvals(R)
            cond = s[0] / s[-1] if s[-1] > 0 else np.inf

            if cond > 1e12:
                return {
                    'coefficients': np.full(k, np.nan),
                    'errors': np.full(k, np.nan),
                    'chi2_red': np.nan,
                    'birge_ratio': np.nan,
                    'condition_number': cond,
                    'status': 'SINGULAR'
                }

            # Solve R β = Qᵀ yw
            qty = Q.T @ yw
            beta = linalg.solve_triangular(R, qty)

            # Residuals and Statistics
            y_pred = X @ beta
            residuals = y - y_pred
            rss = np.sum(weights * residuals**2)
            mse = rss / dof

            # Formal Covariance Matrix: (XᵀWX)⁻¹ = (RᵀR)⁻¹
            # (RᵀR)⁻¹ = R⁻¹ (R⁻ᵀ)
            R_inv = linalg.inv(R)
            cov_raw = (R_inv @ R_inv.T)
            errors_formal = np.sqrt(np.diag(cov_raw))

            # Scale by sqrt(MSE) to get standard errors (CRITICAL FIX)
            # Standard formula: error = sqrt(diag((X'X)^-1)) * sqrt(MSE)
            chi2_red = rss / dof
            mse_unbiased = chi2_red  # Since we're using standardized residuals
            errors_ols = errors_formal * np.sqrt(mse_unbiased)

            # Birge Scaling (only scale up if chi2_red > 1)
            birge_ratio = np.sqrt(chi2_red)
            scaling_factor = max(1.0, birge_ratio) if scale_errors_by_birge else 1.0

            errors_birge = errors_ols * scaling_factor
            cov_scaled = cov_raw * (mse_unbiased * scaling_factor**2)

        return {
            'coefficients': beta,
            'errors': errors_birge,
            'chi2_red': chi2_red,
            'birge_ratio': birge_ratio,
            'condition_number': cond,
            'n_obs': n,
            'dof': dof,
            'rss': rss,
            'mse': mse_unbiased,
            'cov': cov_scaled
        }
    except (linalg.LinAlgError, ValueError):
        return {'coefficients': np.full(k, np.nan), 'errors': np.full(k, np.nan), 
                'chi2_red': np.nan, 'birge_ratio': np.nan, 'condition_number': np.nan}

def linear_regression(y: np.ndarray, x: np.ndarray, weights: np.ndarray = None) -> Dict:
    """
    Backward-compatible wrapper for robust_regression.
    Performs y = A*x + B with Birge error scaling.
    """
    n = len(y)
    X = np.column_stack([x, np.ones(n)])
    res = robust_regression(y, X, weights=weights)
    
    A, B = res['coefficients']
    A_err, B_err = res['errors']
    
    # Convert to η
    eta = A / ETA_SCALE_FACTOR
    eta_err = A_err / ETA_SCALE_FACTOR
    
    # Logging
    from scripts.utils.logger import print_status, get_verbose_mode
    if get_verbose_mode() and np.isfinite(eta):
        print_status(f"Linear Regression Complete (N={n}, DOF={res['dof']}):", "CALC")
        print_status(f"  RSS={res['rss']:.6e}, MSE={res['mse']:.6e}, χ²_red={res['chi2_red']:.3f}", "CALC")
        print_status(f"  Birge Ratio = {res['birge_ratio']:.3f} (Scaling factor applied: {max(1.0, res['birge_ratio']):.3f})", "CALC")
        print_status(f"  Condition Number κ(R) = {res['condition_number']:.2e}", "CALC")
        print_status(f"  Final η = {eta:.8e} ± {eta_err:.8e}", "CALC")
        
    return {
        'amplitude': A,
        'amplitude_error': A_err,
        'intercept': B,
        'intercept_error': B_err,
        'eta': eta,
        'eta_error': eta_err,
        'chi2_red': res['chi2_red'],
        'birge_ratio': res['birge_ratio'],
        'condition_number': res['condition_number'],
        'rss': res['rss'],
        'n_obs': n,
        'method': 'QR-Decomposition with Birge Scaling',
        'regression_metrics': res
    }

def cluster_robust_variance(X: np.ndarray, residuals: np.ndarray,
                          cluster_ids: np.ndarray,
                          small_sample_correction: bool = True) -> Dict:
    """
    Compute cluster-robust (sandwich) covariance matrix and standard errors.

    Uses the Liang-Zeger sandwich estimator with an optional finite-cluster
    correction (G/(G-1)).  Suitable for panel or grouped data where
    observations within a cluster may be correlated.

    Parameters
    ----------
    X : np.ndarray, shape (n, k)
        Design matrix (including intercept column if needed).
    residuals : np.ndarray, shape (n,)
        Residuals from the regression model.
    cluster_ids : np.ndarray, shape (n,)
        Cluster identifier for each observation.
    small_sample_correction : bool, default True
        If True, multiply the meat matrix by G/(G-1) (Cameron-Miller
        finite-cluster correction).

    Returns
    -------
    dict with keys:
        - cov_cluster: cluster-robust covariance matrix, shape (k, k)
        - se_cluster: cluster-robust standard errors, shape (k,)
        - n_clusters: number of distinct clusters
    """
    n, k = X.shape
    XtX = X.T @ X
    XtX_inv = np.linalg.pinv(XtX, rcond=1e-10, hermitian=True)

    unique_clusters = np.unique(cluster_ids)
    G = int(len(unique_clusters))

    meat = np.zeros((k, k))
    for g in unique_clusters:
        mask = cluster_ids == g
        Xg = X[mask]
        ug = residuals[mask]
        # Numerically stable: (Xg * ug[:, None]).T @ (Xg * ug[:, None])
        # avoids overflow from np.outer(ug, ug) with large residuals.
        Xg_ug = Xg * ug[:, np.newaxis]
        meat += Xg_ug.T @ Xg_ug

    cov_cluster = XtX_inv @ meat @ XtX_inv

    if small_sample_correction and G > 1:
        cov_cluster *= G / (G - 1)

    se_cluster = np.sqrt(np.diag(cov_cluster))

    return {
        'cov_cluster': cov_cluster,
        'se_cluster': se_cluster,
        'n_clusters': G
    }


def weighted_linear_regression(residuals: np.ndarray, cos_elong: np.ndarray,
                               weights: np.ndarray, logger=None) -> Dict:
    """
    Deprecated: Use linear_regression or robust_regression instead.
    Maintained for backward compatibility.
    """
    return linear_regression(residuals, cos_elong, weights)


def detect_outliers_iqr(residuals: np.ndarray, multiplier: float = 3.0) -> np.ndarray:
    """
    Detect outliers using IQR (Interquartile Range) method.

    Args:
        residuals: Array of residual values
        multiplier: IQR multiplier for outlier threshold (default 3.0)

    Returns:
        Boolean mask where True indicates outlier
    """
    q25 = np.percentile(residuals, 25)
    q75 = np.percentile(residuals, 75)
    iqr = q75 - q25
    lower_bound = q25 - multiplier * iqr
    upper_bound = q75 + multiplier * iqr

    from scripts.utils.logger import print_status, get_verbose_mode
    if get_verbose_mode():
        mask = (residuals < lower_bound) | (residuals > upper_bound)
        print_status(
            f"IQR Outlier Detection: Q1={q25:.3f}, Q3={q75:.3f}, IQR={iqr:.3f}", "CALC")
        print_status(
            f"  Bounds: [{lower_bound:.3f}, {upper_bound:.3f}], Count={np.sum(mask)}", "CALC")

    return (residuals < lower_bound) | (residuals > upper_bound)


def detect_outliers_sigma(residuals: np.ndarray, sigma_threshold: float = 5.0) -> np.ndarray:
    """
    Detect outliers using Median Absolute Deviation (MAD) method with sigma threshold.
    
    Methodological Justification:
    ----------------------------
    - MAD-based outlier detection is robust against heavy-tailed distributions
    - 6σ threshold corresponds to ~1 in 500 million for Gaussian distribution
    - For 26,207 observations, expected false positives: ~0.05 (essentially zero)
    - This conservative threshold removes only extreme outliers while preserving data integrity
    - Alternative thresholds evaluated: 5σ (too aggressive), 8σ (too permissive)
    - 6σ selected as optimal balance between outlier removal and data preservation
    
    Statistical Rationale:
    ---------------------
    - Converts MAD to standard deviation: σ ≈ 1.4826 × MAD (for normal distribution)
    - Threshold = 6 × 1.4826 × MAD ≈ 8.9 × MAD
    - This removes observations beyond ~8.9 median absolute deviations
    - Equivalent to ~6 standard deviations for Gaussian data
    
    References:
    -----------
    - Leys et al. 2013, "Identifying Outliers: Do Not Use Standard Deviation"
    - Rousseeuw & Leroy 1987, "Robust Regression and Outlier Detection"
    - Hoaglin et al. 1983, "Understanding Robust and Exploratory Data Analysis"

    Args:
        residuals: Array of residual values [meters]
        sigma_threshold: Number of standard deviations for outlier threshold (default 5.0)
                        Pipeline standard: 6.0 (used in steps 002, 003, 005, 016, 024, 025, 052)

    Returns:
        Boolean mask where True indicates outlier (should be removed)
    """
    median = np.median(residuals)
    mad = np.median(np.abs(residuals - median))
    # Convert MAD to standard deviation (for normal distribution, sigma ≈ 1.4826 * MAD)
    sigma = 1.4826 * mad
    threshold = sigma_threshold * sigma

    from scripts.utils.logger import print_status, get_verbose_mode
    if get_verbose_mode():
        mask = np.abs(residuals - median) > threshold
        print_status(
            f"Sigma Outlier Detection: median={median:.3f}, MAD={mad:.3f}, sigma={sigma:.3f}", "CALC")
        print_status(
            f"  Threshold={threshold:.3f} ({sigma_threshold}σ), Count={np.sum(mask)}", "CALC")

    return np.abs(residuals - median) > threshold


def detect_outliers_isolation_forest(residuals: np.ndarray, elongation: np.ndarray,
                                     contamination: float = 0.01) -> np.ndarray:
    """
    Detect outliers using Isolation Forest (requires sklearn).
    Falls back to sigma method if sklearn not available.

    Args:
        residuals: Array of residual values
        elongation: Array of elongation values (for multivariate analysis)
        contamination: Expected proportion of outliers (default 0.01)

    Returns:
        Boolean mask where True indicates outlier
    """
    try:
        from sklearn.ensemble import IsolationForest

        # Stack residuals and elongation for multivariate analysis
        X = np.column_stack([residuals, elongation])

        clf = IsolationForest(contamination=contamination,
                              random_state=TEP_CONFIG.get("RANDOM_SEED", 42), n_jobs=-1)
        outlier_pred = clf.fit_predict(X)
        return outlier_pred == -1  # -1 indicates outlier
    except ImportError:
        raise ImportError(
            "scikit-learn is required for Isolation Forest outlier detection. "
            "Install it via 'pip install scikit-learn' or choose a different "
            "methodology explicitly. Silent fallbacks are prohibited."
        )


def compute_minimum_detectable_eta(
    n_obs: int,
    residual_rms: float,
    sigma_cos_elong: float,
    confidence_level: float = 3.0,
) -> dict:
    """
    Compute the minimum detectable Nordtvedt parameter.

    For OLS of residuals y on x = cos(elongation), Pearson r and slope A satisfy
    r = A * sigma_x / sigma_y (population analogue). Sensitivity must use both
    sigma_y (residual scale) and sigma_x (predictor std).

    Args:
        n_obs: Number of observations
        residual_rms: RMS (or typical scale) of residuals in meters [sigma_y]
        sigma_cos_elong: Standard deviation of cos(elongation) over the sample [sigma_x]
        confidence_level: Confidence level in sigma (default 3.0)

    Returns:
        Dictionary with sensitivity metrics
    """
    if n_obs <= 0:
        raise ValueError(f"n_obs must be positive, got {n_obs}")
    if residual_rms <= 0:
        raise ValueError(f"residual_rms must be positive, got {residual_rms}")
    if sigma_cos_elong <= 0:
        raise ValueError(
            f"sigma_cos_elong must be positive, got {sigma_cos_elong}; "
            "pass e.g. float(np.std(cos_elong)) from the same dataset."
        )

    # Standard error of correlation coefficient (null, large-n)
    sigma_r = 1.0 / np.sqrt(n_obs)

    # Minimum detectable correlation at given confidence level
    r_min = confidence_level * sigma_r

    # A = r * sigma_y / sigma_x  =>  A_min = r_min * residual_rms / sigma_cos_elong
    A_min = r_min * residual_rms / sigma_cos_elong

    # Convert to minimum detectable eta (A = 13 * eta)
    eta_min = A_min / ETA_SCALE_FACTOR

    # Power analysis for specific eta values
    eta_test_values = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
    power_results = {}

    for eta_test in eta_test_values:
        A_test = ETA_SCALE_FACTOR * eta_test
        r_test = A_test * sigma_cos_elong / residual_rms
        z_score = r_test / sigma_r
        # Two-tailed power at α=0.05 (Z_α/2 = 1.96)
        power = 1 - (stats.norm.cdf(Z_ALPHA_2 - z_score) -
                     stats.norm.cdf(-Z_ALPHA_2 - z_score))
        # Clamp to [0,1] to handle numerical edge cases
        power = max(0.0, min(1.0, power))
        power_results[f'eta_{eta_test:.0e}'] = {
            'eta': eta_test,
            'amplitude_m': A_test,
            'correlation_r': r_test,
            'z_score': z_score,
            'power': power
        }

    return {
        'n_observations': n_obs,
        'residual_rms_m': residual_rms,
        'sigma_cos_elong': float(sigma_cos_elong),
        'confidence_level_sigma': confidence_level,
        'sigma_r': sigma_r,
        'r_min': r_min,
        'A_min_m': A_min,
        'eta_min': eta_min,
        'power_analysis': power_results
    }


def compute_sensitivity_by_sample_size(
    residual_rms: float,
    sigma_cos_elong: float,
    eta_target: float = 1e-4,
    sample_sizes: list = None,
) -> dict:
    """
    Compute required sample size to detect a given eta with 95% power.

    Args:
        residual_rms: RMS of residuals in meters [sigma_y]
        sigma_cos_elong: Std. dev. of cos(elongation) [sigma_x]
        eta_target: Target Nordtvedt parameter to detect
        sample_sizes: List of sample sizes to test (default: log-spaced from 100 to 100000)

    Returns:
        Dictionary with sensitivity vs sample size
    """
    if residual_rms <= 0:
        raise ValueError(f"residual_rms must be positive, got {residual_rms}")
    if sigma_cos_elong <= 0:
        raise ValueError(f"sigma_cos_elong must be positive, got {sigma_cos_elong}")
    if eta_target <= 0:
        raise ValueError(f"eta_target must be positive, got {eta_target}")

    if sample_sizes is None:
        sample_sizes = np.logspace(2, 5, 20).astype(int)

    results = {}

    for n in sample_sizes:
        sigma_r = 1.0 / np.sqrt(n)
        A_target = ETA_SCALE_FACTOR * eta_target
        r_target = A_target * sigma_cos_elong / residual_rms
        z_score = r_target / sigma_r
        # Two-tailed power at α=0.05 (Z_α/2 = 1.96)
        power = 1 - (stats.norm.cdf(Z_ALPHA_2 - z_score) -
                     stats.norm.cdf(-Z_ALPHA_2 - z_score))
        # Clamp to [0,1] to handle numerical edge cases
        power = max(0.0, min(1.0, power))

        results[int(n)] = {
            'n': int(n),
            'sigma_r': sigma_r,
            'z_score': z_score,
            'power': power
        }

    # Find minimum sample size for 95% power
    n_for_95_power = None
    for n in sorted(results.keys()):
        if results[n]['power'] >= 0.95:
            n_for_95_power = n
            break

    return {
        'eta_target': eta_target,
        'residual_rms_m': residual_rms,
        'sigma_cos_elong': float(sigma_cos_elong),
        'sensitivity_by_n': results,
        'n_for_95_percent_power': n_for_95_power
    }


def compute_sensitivity_by_precision(
    n_obs: int,
    sigma_cos_elong: float,
    eta_target: float = 1e-4,
    rms_values: list = None,
) -> dict:
    """
    Compute required precision to detect a given eta with 95% power.

    Args:
        n_obs: Number of observations
        sigma_cos_elong: Std. dev. of cos(elongation) [sigma_x]
        eta_target: Target Nordtvedt parameter to detect
        rms_values: List of RMS values to test in meters (default: 0.01 to 0.5 m)

    Returns:
        Dictionary with sensitivity vs precision
    """
    if n_obs <= 0:
        raise ValueError(f"n_obs must be positive, got {n_obs}")
    if sigma_cos_elong <= 0:
        raise ValueError(f"sigma_cos_elong must be positive, got {sigma_cos_elong}")
    if eta_target <= 0:
        raise ValueError(f"eta_target must be positive, got {eta_target}")

    if rms_values is None:
        rms_values = np.linspace(0.01, 0.5, 20)

    results = {}

    for rms in rms_values:
        sigma_r = 1.0 / np.sqrt(n_obs)
        A_target = ETA_SCALE_FACTOR * eta_target
        r_target = A_target * sigma_cos_elong / rms
        z_score = r_target / sigma_r
        # Two-tailed power at α=0.05 (Z_α/2 = 1.96)
        power = 1 - (stats.norm.cdf(Z_ALPHA_2 - z_score) -
                     stats.norm.cdf(-Z_ALPHA_2 - z_score))
        # Clamp to [0,1] to handle numerical edge cases
        power = max(0.0, min(1.0, power))

        results[float(rms)] = {
            'rms_m': rms,
            'sigma_r': sigma_r,
            'z_score': z_score,
            'power': power
        }

    # Find maximum RMS for 95% power
    rms_for_95_power = None
    for rms in sorted(results.keys(), reverse=True):
        if results[rms]['power'] >= 0.95:
            rms_for_95_power = rms
            break

    return {
        'eta_target': eta_target,
        'n_observations': n_obs,
        'sigma_cos_elong': float(sigma_cos_elong),
        'sensitivity_by_rms': results,
        'rms_for_95_percent_power_m': rms_for_95_power
    }


def fishers_combined_probability(p_values: np.ndarray) -> Dict:
    """
    Fisher's combined probability test for aggregating independent p-values.

    Fisher's method combines p-values from multiple independent statistical
    tests into a single test statistic. This is useful for meta-analysis when
    combining results from different stations or analysis methods.

    The test statistic X² = -2 Σ ln(p_i) follows a chi-squared distribution
    with 2k degrees of freedom, where k is the number of tests.

    Args:
        p_values: Array of p-values from independent tests (must be in (0, 1])

    Returns:
        Dictionary with:
        - combined_chi2: Fisher's combined test statistic
        - df: Degrees of freedom (2 * number of p-values)
        - combined_p: Combined p-value
        - n_tests: Number of tests combined
        - method: 'Fisher's combined probability test'

    References:
        Fisher RA (1948) Combining independent tests of significance.
        American Statistician 2:30.
    """
    p_values = np.asarray(p_values, dtype=float)
    n_tests = len(p_values)

    if n_tests == 0:
        return {
            'combined_chi2': np.nan,
            'df': 0,
            'combined_p': np.nan,
            'n_tests': 0,
            'method': "Fisher's combined probability test"
        }

    # Check for invalid p-values
    if np.any(p_values <= 0) or np.any(p_values > 1):
        raise ValueError("p-values must be in (0, 1]")

    # Fisher's method: X² = -2 Σ ln(p_i) ~ χ²(2k)
    # Use log transform to avoid underflow
    log_p_sum = np.sum(np.log(p_values))
    chi2_stat = -2.0 * log_p_sum
    df = 2 * n_tests

    # Combined p-value from chi-squared distribution
    combined_p = 1 - stats.chi2.cdf(chi2_stat, df)

    return {
        'combined_chi2': float(chi2_stat),
        'df': int(df),
        'combined_p': float(combined_p),
        'n_tests': int(n_tests),
        'method': "Fisher's combined probability test"
    }


def steiger_z_test(r1: float, r2: float, n: int, r12: float = None) -> Dict:
    """
    Steiger's Z-test for comparing two correlated correlation coefficients.

    When comparing correlations from the same sample (e.g., correlation of
    residuals with cos(D) vs correlation with some systematic variable),
    the correlations are not independent. Steiger's test accounts for this
    correlation between correlations.

    This is distinct from Fisher's Z-test for independent correlations.

    Args:
        r1: First correlation coefficient
        r2: Second correlation coefficient
        n: Sample size
        r12: Correlation between the two variables being correlated with
             the common variable (default: 0 for orthogonal predictors)

    Returns:
        Dictionary with:
        - z_stat: Steiger's Z statistic
        - p_value: Two-tailed p-value
        - r1, r2: Input correlations
        - r12: Correlation between predictors
        - n: Sample size
        - method: "Steiger's Z-test for dependent correlations"

    References:
        Steiger JH (1980) Tests for comparing elements of a correlation matrix.
        Psychological Bulletin 87:245-251.
    """
    if r12 is None:
        r12 = 0.0  # Assume orthogonal predictors if not specified

    # Fisher Z-transform
    z1 = 0.5 * np.log((1 + r1) / (1 - r1))
    z2 = 0.5 * np.log((1 + r2) / (1 - r2))

    # Mean correlation for covariance calculation
    r_mean = (r1 + r2) / 2.0

    # Covariance between correlations (Steiger 1980, Eq. 6)
    # cov(r1, r2) = [(r12 * (1 - r1^2 - r2^2 - r12^2)/2 + r1*r2*(1-r1^2-r2^2-r12^2)/2
    #                + r1*r2*r12^2] / (1 - r_mean^2)^2
    # Simplified form assuming orthogonal predictors (r12 = 0):
    # cov(r1, r2) = r1 * r2 * (1 - r1^2 - r2^2) / (2 * (1 - r_mean^2)^2)

    if abs(r_mean) >= 1.0:
        # Degenerate case
        return {
            'z_stat': np.nan,
            'p_value': np.nan,
            'r1': r1,
            'r2': r2,
            'r12': r12,
            'n': n,
            'method': "Steiger's Z-test for dependent correlations (degenerate)"
        }

    # Variance of Fisher z then delta-method to r: Var(r) ≈ (1-r^2)^2 / (n-3)
    denom = max(n - 3, 1)
    var_r1 = (1 - r1**2) ** 2 / denom
    var_r2 = (1 - r2**2) ** 2 / denom

    # Covariance between correlations (Steiger 1980); leading scale ~ 1/n
    cov_r1_r2 = ((r12 * (1 - r1**2 - r2**2 - r12**2)) / 2.0 +
                 r1 * r2 * (1 - r1**2 - r2**2 - r12**2) / 2.0 +
                 r1 * r2 * r12**2) / (1 - r_mean**2)**2 / n

    # Standard error of difference
    se_diff = np.sqrt(var_r1 + var_r2 - 2 * cov_r1_r2)

    if se_diff == 0:
        return {
            'z_stat': np.nan,
            'p_value': np.nan,
            'r1': r1,
            'r2': r2,
            'r12': r12,
            'n': n,
            'method': "Steiger's Z-test for dependent correlations (zero SE)"
        }

    # Z-statistic
    z_stat = (z1 - z2) / se_diff

    # Two-tailed p-value
    p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))

    return {
        'z_stat': float(z_stat),
        'p_value': float(p_value),
        'r1': r1,
        'r2': r2,
        'r12': r12,
        'n': n,
        'se_diff': float(se_diff),
        'method': "Steiger's Z-test for dependent correlations"
    }


def require_step003_eta_ols(payload: dict) -> float:
    """Return finite eta_ols from step_003_statistical_analysis.json (strict)."""
    if "eta_ols" not in payload:
        raise KeyError(
            "step_003_statistical_analysis.json missing required key 'eta_ols'"
        )
    v = payload["eta_ols"]
    if v is None:
        raise ValueError("eta_ols is null")
    out = float(v)
    if not np.isfinite(out):
        raise ValueError("eta_ols must be finite")
    return out


def require_step003_eta_ols_error(payload: dict) -> float:
    """Return positive finite eta_ols_error from step_003 output (strict)."""
    if "eta_ols_error" not in payload:
        raise KeyError(
            "step_003_statistical_analysis.json missing required key 'eta_ols_error'"
        )
    v = payload["eta_ols_error"]
    if v is None:
        raise ValueError("eta_ols_error is null")
    out = float(v)
    if not np.isfinite(out) or out <= 0:
        raise ValueError("eta_ols_error must be finite and positive")
    return out
