"""
Astronomical utility functions for TEP-LLR analysis.

This module contains shared astronomical calculations including elongation,
phase, and related computations that are used across multiple parsing modules.
"""

import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CANONICAL_EPHEMERIS = PROJECT_ROOT / "data" / "raw" / "de440.bsp"


def resolve_skyfield_kernel_path(project_root: Path | None = None) -> Path:
    """Return the manifest-checked JPL kernel path (de440.bsp)."""
    root = project_root or PROJECT_ROOT
    eph_path = root / "data" / "raw" / "de440.bsp"
    if not eph_path.exists():
        raise FileNotFoundError(
            f"Required Skyfield kernel missing: {eph_path}. "
            "Place de440.bsp under data/raw and verify via Step 000 manifest."
        )
    return eph_path


def load_skyfield_planets(project_root: Path | None = None):
    """Load Skyfield planets from the canonical de440 kernel."""
    from skyfield.api import load

    eph_path = resolve_skyfield_kernel_path(project_root)
    return load(str(eph_path)), eph_path


def compute_elongation(df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
    """
    Compute lunar elongation (angular separation Sun–Moon as seen from Earth) using Skyfield/de440.

    Elongation is the geometric separation angle; ``cos(elongation)`` matches the
    sky-plane dot product used with range residuals.

    Args:
        df: DataFrame containing 'date_julian' column
        verbose: If True, print debug information

    Returns:
        DataFrame with added 'elongation_rad' column
    """
    if verbose:
        print_status("Computing precise elongation via Skyfield (de440.bsp)...", "INFO")

    from skyfield.api import load as skyfield_load

    planets, eph_path = load_skyfield_planets()
    if verbose:
        print_status(f"Using ephemeris: {eph_path.relative_to(PROJECT_ROOT)}", "INFO")
    earth = planets['earth']
    moon = planets['moon']
    sun = planets['sun']
    ts = skyfield_load.timescale()

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

    # dot = cos(theta) for true geometric Sun–Moon separation as seen from Earth
    cos_d = np.einsum('ij,ij->j', m_unit, s_unit)
    cos_d = np.clip(cos_d, -1.0, 1.0)

    # Use geometric elongation so cos(elongation) matches the sky-plane cosine;
    # ecliptic longitude difference alone is not the separation angle when latitudes differ.
    elongation_rad = np.arccos(cos_d)

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
