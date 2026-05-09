"""
Astronomical utility functions for TEP-LLR analysis.

This module contains shared astronomical calculations including elongation,
phase, and related computations that are used across multiple parsing modules.
"""

import numpy as np
import pandas as pd

# Import constants from llr_constants module
from scripts.utils.llr_constants import (
    SYNODIC_PERIOD_DAYS,
    REFERENCE_NEW_MOON_JD,
    SYNODIC_PERIOD_RADIANS
)


def compute_elongation(df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
    """
    Compute precise lunar elongation (angular separation from Sun) using Skyfield/DE421.

    Args:
        df: DataFrame containing 'date_julian' column
        verbose: If True, print debug information

    Returns:
        DataFrame with added 'elongation_rad' column
    """
    from skyfield.api import load
    from pathlib import Path
    
    PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
    eph_path = PROJECT_ROOT / "de421.bsp"
    
    if eph_path.exists():
        if verbose:
            print_status(f"Computing precise elongation via Skyfield ({eph_path.name})...", "INFO")
        
        planets = load(str(eph_path))
        earth = planets['earth']
        moon = planets['moon']
        sun = planets['sun']
        ts = load.timescale()
        
        times = ts.tt(jd=df['date_julian'].values)
        
        # Get positions in GCRS
        e_pos = earth.at(times)
        m_pos = e_pos.observe(moon).apparent()
        s_pos = e_pos.observe(sun).apparent()
        
        # Angle between Moon and Sun as seen from Earth
        # Using dot product of unit vectors for stability
        m_vec = m_pos.position.au
        s_vec = s_pos.position.au
        
        m_unit = m_vec / np.linalg.norm(m_vec, axis=0)
        s_unit = s_vec / np.linalg.norm(s_vec, axis=0)
        
        # dot = cos(theta)
        cos_d = np.einsum('ij,ij->j', m_unit, s_unit)
        # Ensure range [-1, 1]
        cos_d = np.clip(cos_d, -1.0, 1.0)
        
        # elongation_rad = arccos(cos_d)
        # Note: arccos gives [0, pi]. We need [0, 2pi] to distinguish waxing/waning.
        # We need to check if Moon is ahead or behind Sun in longitude.
        # Simplification: TEP signal depends on cos(D), so arccos is technically sufficient
        # if we only care about the cos(D) predictor. 
        # But for completeness, we use the phase-angle logic.
        
        # Get ecliptic longitudes
        m_lon = m_pos.ecliptic_latlon()[1].radians
        s_lon = s_pos.ecliptic_latlon()[1].radians
        
        elongation_rad = (m_lon - s_lon) % (2 * np.pi)
    else:
        print_status("Warning: de421.bsp not found. Falling back to mean synodic approximation.", "WARNING")
        # Calculate phase in cycles based on synodic period
        phase_cycles = (df['date_julian'] - REFERENCE_NEW_MOON_JD) / SYNODIC_PERIOD_DAYS
        elongation_rad = (phase_cycles * SYNODIC_PERIOD_RADIANS) % SYNODIC_PERIOD_RADIANS

    df['elongation_rad'] = elongation_rad

    if verbose:
        print_status(f"Elongation computation complete (N={len(df)}):", "CALC")
        print_status(f"  Range: {elongation_rad.min():.6f} to {elongation_rad.max():.6f} rad", "CALC")

    return df


def compute_cos_elongation(df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
    """
    Compute cosine of elongation for TEP signal analysis.

    The TEP signal is proportional to cos(elongation), so this function
    computes the cosine of the elongation angle.

    Args:
        df: DataFrame containing 'elongation_rad' column
        verbose: If True, print debug information

    Returns:
        DataFrame with added 'cos_elong' column containing cos(elongation)
    """
    if 'elongation_rad' not in df.columns:
        raise ValueError("DataFrame must contain 'elongation_rad' column")

    df['cos_elong'] = np.cos(df['elongation_rad'])

    from scripts.utils.logger import get_verbose_mode
    if get_verbose_mode():
        print_status("Cos(elongation) computation complete:", "CALC")
        print_status(
            f"  Range: {df['cos_elong'].min():.6f} to {df['cos_elong'].max():.6f}", "CALC")
        print_status(
            f"  Mean:  {df['cos_elong'].mean():.6f}, RMS: {np.sqrt(np.mean(df['cos_elong']**2)):.6f}", "CALC")

    return df


def get_phase_from_elongation(elongation_rad: float) -> str:
    """
    Get lunar phase description from elongation angle.

    Args:
        elongation_rad: Elongation angle in radians

    Returns:
        String description of lunar phase
    """
    # Normalize to [0, 2π)
    elongation_rad = elongation_rad % (2 * np.pi)

    if elongation_rad < 0.1 or elongation_rad > 2 * np.pi - 0.1:
        return "New Moon"
    elif elongation_rad < np.pi - 0.1:
        return "Waxing Crescent" if elongation_rad < np.pi / 2 else "Waxing Gibbous"
    elif elongation_rad > np.pi + 0.1 and elongation_rad < 2 * np.pi - 0.1:
        return "Waning Gibbous" if elongation_rad < 3 * np.pi / 2 else "Waning Crescent"
    else:
        return "Full Moon"


def mask_phase_region(elongation_rad: np.ndarray,
                      center_rad: float,
                      width_rad: float = 0.5) -> np.ndarray:
    """
    Create a boolean mask for observations within a phase region.

    Args:
        elongation_rad: Array of elongation angles in radians
        center_rad: Center of the phase region in radians
        width_rad: Half-width of the phase region in radians (default: 0.5 rad ≈ 28.6°)

    Returns:
        Boolean mask for observations within the phase region
    """
    # Normalize elongation to [0, 2π)
    elongation_normalized = elongation_rad % (2 * np.pi)
    center_normalized = center_rad % (2 * np.pi)

    # Handle wrap-around case (e.g., near 0/2π boundary)
    if center_normalized < width_rad:
        # Region wraps around 0/2π
        mask = (elongation_normalized < center_normalized + width_rad) | \
               (elongation_normalized > 2 * np.pi - (width_rad - center_normalized))
    elif center_normalized > 2 * np.pi - width_rad:
        # Region wraps around 2π/0
        mask = (elongation_normalized > center_normalized - width_rad) | \
               (elongation_normalized < width_rad - (2 * np.pi - center_normalized))
    else:
        # Normal case
        mask = np.abs(elongation_normalized - center_normalized) < width_rad

    return mask


# Import print_status from logger module
from scripts.utils.logger import print_status
