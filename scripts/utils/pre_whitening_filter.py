#!/usr/bin/env python3
"""
Pre-Whitening Filter for TEP-LLR
Removes dominant non-synodic seasonal and systematic harmonics from residuals.
"""
import numpy as np
import pandas as pd
from typing import List
from scripts.utils.llr_constants import SYNODIC_PERIOD_DAYS
from scripts.utils.logger import print_status


def identify_peak_harmonics(phase: np.ndarray, y: np.ndarray,
                            n_peaks: int = 5,
                            exclude_synodic_factor: float = 0.05) -> List[float]:
    """
    Identifies the top N frequencies with highest power, excluding the synodic region.

    Parameters
    ----------
    phase : np.ndarray
        The phase angle in radians (e.g., elongation_rad). Must match the phase
        used in the downstream analysis so that pre-whitening and detection share
        the same basis. Using uniform time (e.g. date_julian) will cause a basis
        mismatch because lunar elongation is not uniform in time.
    y : np.ndarray
        Residual values [m].

    WARNING: This is a data-driven peak picker. If the true signal has sidebands
    near the synodic frequency, this filter may inadvertently remove or attenuate
    genuine TEP signal power. Use with caution; verify null-test behaviour when
    applying pre-whitening to detection steps.
    """
    # Use a coarse frequency sweep (factors of synodic frequency)
    # Range: 0.1 to 4.0 * synodic
    factors = np.linspace(0.1, 4.0, 400)

    snrs = []
    for f in factors:
        # Exclude synodic region to avoid over-whitening the signal we want to detect
        if abs(f - 1.0) < exclude_synodic_factor:
            snrs.append(0)
            continue

        # Fit sin/cos pair for this frequency in the phase domain
        # phase is already in radians, so basis is cos(f * phase), sin(f * phase)
        cos_term = np.cos(f * phase)
        sin_term = np.sin(f * phase)

        X = np.column_stack([cos_term, sin_term, np.ones_like(phase)])
        try:
            coeffs, res_sum, rank, s = np.linalg.lstsq(X, y, rcond=None)
            rss = res_sum[0] if len(
                res_sum) > 0 else np.sum((y - X @ coeffs)**2)
            n = len(y)
            sigma2 = rss / (n - 3)
            cov = sigma2 * np.linalg.inv(X.T @ X)
            # Amplitude squared A^2 = c1^2 + s1^2
            # For simplicity, we use the joint significance
            snr = np.sqrt(coeffs[0]**2 + coeffs[1]**2) / \
                np.sqrt(cov[0, 0] + cov[1, 1])
            snrs.append(snr)
        except (np.linalg.LinAlgError, ValueError):
            snrs.append(0)

    # Find peaks
    peaks_idx = np.argsort(snrs)[-n_peaks:][::-1]
    peak_factors = [factors[i] for i in peaks_idx]

    from scripts.utils.logger import get_verbose_mode
    if get_verbose_mode():
        print_status(
            f"Spectral analysis complete (N_scan={len(factors)}):", "CALC")
        for i, idx in enumerate(peaks_idx):
            print_status(
                f"  Peak {i+1}: f={factors[idx]:.4f} * synodic, SNR={snrs[idx]:.2f}", "CALC")

    return peak_factors


def apply_pre_whitening(df: pd.DataFrame, n_harmonics: int = 5, verbose: bool = False,
                        extra_harmonics: List[float] = None) -> pd.DataFrame:
    """
    Subtracts the top N non-synodic harmonics from the residuals.
    Returns a copy of the dataframe with 'residual_whitened_m' column.
    """
    df_clean = df.copy()
    y = df_clean['residual_m'].values

    # Use elongation_rad as the phase basis so pre-whitening shares the same
    # basis as the downstream analysis (null test, primary regression, etc.).
    # Using uniform time (date_julian) causes a basis mismatch because lunar
    # elongation is not uniform in time (orbital eccentricity + perturbations).
    if 'elongation_rad' in df_clean.columns:
        phase = df_clean['elongation_rad'].values
    else:
        print_status(
            "WARNING: elongation_rad not found in dataframe. Falling back to date_julian "
            "for pre-whitening basis. This may cause basis mismatch with downstream analysis.",
            "WARNING")
        t = df_clean['date_julian'].values
        phase = (t - np.min(t)) * (2 * np.pi / SYNODIC_PERIOD_DAYS)

    peaks = identify_peak_harmonics(phase, y, n_peaks=n_harmonics)

    if extra_harmonics:
        for f in extra_harmonics:
            if not any(abs(f - p) < 1e-6 for p in peaks):
                peaks.append(float(f))

    if verbose:
        print_status(
            f"Found {len(peaks)} non-physical harmonics for pre-whitening:", "PROCESS")
        for p in peaks:
            print_status(
                f"  Factor {p:.4f} * Synodic (~{SYNODIC_PERIOD_DAYS/p:.1f} days)", "INFO")

    y_whitened = y.copy()

    # Construct the full design matrix for all harmonics simultaneously
    # Column 0 is the constant (intercept)
    X_full = [np.ones_like(phase)]

    for f in peaks:
        X_full.append(np.cos(f * phase))
        X_full.append(np.sin(f * phase))

    X_full = np.column_stack(X_full)

    # Solve for all harmonics jointly
    try:
        coeffs, _, _, _ = np.linalg.lstsq(X_full, y, rcond=None)

        # The first coefficient is the mean/intercept
        # The others are [cos_f1, sin_f1, cos_f2, sin_f2, ...]

        # Calculate the total fit
        y_fit = X_full @ coeffs

        # The whitened residuals are the original minus the joint fit
        y_whitened = y - y_fit

        if verbose:
            print_status(f"Jointly whitened {len(peaks)} harmonics (N_params={X_full.shape[1]})", "SUCCESS")
    except np.linalg.LinAlgError as e:
        print_status(f"Error in joint whitening fit: {e}. Falling back to iterative.", "WARNING")
        for f in peaks:
            c = np.cos(f * phase)
            s = np.sin(f * phase)
            M = np.column_stack([c, s, np.ones_like(phase)])
            cf, _, _, _ = np.linalg.lstsq(M, y_whitened, rcond=None)
            y_whitened -= (cf[0] * c + cf[1] * s)

    df_clean['residual_whitened_m'] = y_whitened
    return df_clean


def validate_pre_whitening(df: pd.DataFrame, verbose: bool = False) -> dict:
    """
    Validate that pre-whitening does not artificially enhance or suppress
    the target synodic signal.

    Computes Pearson r(residuals, cos(elongation)) before and after whitening.
    A genuine signal should either persist or weaken slightly after removing
    nuisance harmonics; a large increase would indicate the filter is creating
    artificial orthogonality.
    """
    from scipy import stats

    cos_elong = np.cos(df['elongation_rad'].values)
    y_raw = df['residual_m'].values

    r_raw, p_raw = stats.pearsonr(y_raw, cos_elong)

    df_white = apply_pre_whitening(df, n_harmonics=5, verbose=False)
    y_white = df_white['residual_whitened_m'].values
    r_white, p_white = stats.pearsonr(y_white, cos_elong)

    delta_r = r_white - r_raw
    delta_r_percent = (delta_r / r_raw * 100) if r_raw != 0 else 0.0

    # Flag if whitening changes correlation magnitude by more than 20%
    # or flips the sign unexpectedly
    sign_changed = (r_raw < 0 and r_white > 0) or (r_raw > 0 and r_white < 0)
    magnitude_shift = abs(delta_r_percent)
    flagged = sign_changed or magnitude_shift > 20

    if verbose:
        print_status("PRE-WHITENING VALIDATION", "PROCESS")
        print_status(f"  Raw  r = {r_raw:.6e} (p={p_raw:.2e})", "CALC")
        print_status(f"  White r = {r_white:.6e} (p={p_white:.2e})", "CALC")
        print_status(f"  Δr = {delta_r:.6e} ({delta_r_percent:+.1f}%)", "CALC")
        print_status(f"  Flagged: {flagged}", "WARNING" if flagged else "SUCCESS")

    return {
        'r_raw': float(r_raw),
        'p_raw': float(p_raw),
        'r_whitened': float(r_white),
        'p_whitened': float(p_white),
        'delta_r': float(delta_r),
        'delta_r_percent': float(delta_r_percent),
        'sign_changed': bool(sign_changed),
        'magnitude_shift_percent': float(magnitude_shift),
        'flagged': bool(flagged),
        'interpretation': (
            'WARNING: Pre-whitening significantly altered synodic correlation. '
            'Filter may be overfitting or removing genuine signal components.'
            if flagged else
            'PASS: Synodic correlation stable after pre-whitening.'
        )
    }


if __name__ == "__main__":
    # Test script
    from scripts.utils.logger import TEPLogger
    logger = TEPLogger("pre_whitening_test", verbose=True)

    # Mock data
    t = np.linspace(0, 1000, 2000)
    # TEP signal at f=1.0
    y = 0.01 * np.cos(2 * np.pi * 1.0 / 29.53 * t)
    # Noise artifact at f=1.23
    y += 0.05 * np.cos(2 * np.pi * 1.23 / 29.53 * t + 0.5)
    # White noise
    y += np.random.normal(0, 0.1, len(t))

    # Compute mock elongation_rad from uniform time (only for test data)
    phase = (t - np.min(t)) * (2 * np.pi / 29.53)
    df = pd.DataFrame({'date_julian': t, 'residual_m': y, 'elongation_rad': phase})
    df_whitened = apply_pre_whitening(df, n_harmonics=1, verbose=True)

    print_status(f"Original RMS: {np.std(y):.4f}")
    print_status(
        f"Whitened RMS: {np.std(df_whitened['residual_whitened_m']):.4f}")
