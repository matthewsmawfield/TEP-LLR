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


def identify_peak_harmonics(t: np.ndarray, y: np.ndarray,
                            n_peaks: int = 5,
                            exclude_synodic_factor: float = 0.05) -> List[float]:
    """
    Identifies the top N frequencies with highest power, excluding the synodic region.
    """
    # Use a coarse frequency sweep (factors of synodic frequency)
    # Range: 0.1 to 10.0 * synodic
    factors = np.linspace(0.1, 4.0, 400)
    synodic_freq = 1.0 / SYNODIC_PERIOD_DAYS

    snrs = []
    for f in factors:
        # Exclude synodic region to avoid over-whitening the signal we want to detect
        if abs(f - 1.0) < exclude_synodic_factor:
            snrs.append(0)
            continue

        # Fit sin/cos pair for this frequency
        omega = 2 * np.pi * f * synodic_freq
        cos_term = np.cos(omega * t)
        sin_term = np.sin(omega * t)

        X = np.column_stack([cos_term, sin_term, np.ones_like(t)])
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
    t = df_clean['date_julian'].values
    y = df_clean['residual_m'].values

    # Normalize time to start at 0 for numeric stability
    t_norm = t - np.min(t)

    peaks = identify_peak_harmonics(t_norm, y, n_peaks=n_harmonics)

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
    synodic_freq = 1.0 / SYNODIC_PERIOD_DAYS

    # Construct the full design matrix for all harmonics simultaneously
    # Column 0 is the constant (intercept)
    X_full = [np.ones_like(t_norm)]
    
    for f in peaks:
        omega = 2 * np.pi * f * synodic_freq
        X_full.append(np.cos(omega * t_norm))
        X_full.append(np.sin(omega * t_norm))
        
    X_full = np.column_stack(X_full)
    
    # Solve for all harmonics jointly
    try:
        coeffs, _, _, _ = np.linalg.lstsq(X_full, y, rcond=None)
        
        # The first coefficient is the mean/intercept
        # The others are [cos_f1, sin_f1, cos_f2, sin_f2, ...]
        
        # Calculate the total fit (excluding the intercept if we want to keep the mean)
        # Or better: keep the residuals centered. 
        y_fit = X_full @ coeffs
        
        # The whitened residuals are the original minus the joint fit
        y_whitened = y - y_fit
        
        if verbose:
            print_status(f"Jointly whitened {len(peaks)} harmonics (N_params={X_full.shape[1]})", "SUCCESS")
    except np.linalg.LinAlgError as e:
        print_status(f"Error in joint whitening fit: {e}. Falling back to iterative.", "WARNING")
        # (Fall back logic omitted for brevity, but joint is standard for this N)
        for f in peaks:
            omega = 2 * np.pi * f * synodic_freq
            c = np.cos(omega * t_norm)
            s = np.sin(omega * t_norm)
            M = np.column_stack([c, s, np.ones_like(t_norm)])
            cf, _, _, _ = np.linalg.lstsq(M, y_whitened, rcond=None)
            y_whitened -= (cf[0] * c + cf[1] * s)

    df_clean['residual_whitened_m'] = y_whitened
    return df_clean


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

    df = pd.DataFrame({'date_julian': t, 'residual_m': y})
    df_whitened = apply_pre_whitening(df, n_harmonics=1, verbose=True)

    print_status(f"Original RMS: {np.std(y):.4f}")
    print_status(
        f"Whitened RMS: {np.std(df_whitened['residual_whitened_m']):.4f}")
