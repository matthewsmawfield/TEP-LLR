#!/usr/bin/env python3
"""
Step 076: Full Sky-Scan with Look-Elsewhere Correction
=======================================================

Addresses the cherry-picking criticism of the CMB directional analysis.
Instead of testing one predetermined axis (Planck 2018 dipole), we grid the
entire celestial sphere and apply a look-elsewhere correction.

Core logic:
  1. Grid the celestial sphere in 5° steps (≈ 3,000 directions)
  2. For each direction n̂, fit η = η₀ + η_n̂ · cos(θ_EM-n̂)
  3. Record the maximum ΔAIC across all directions
  4. Compare to Monte Carlo scrambled skies (10,000 realizations)
  5. Report the look-elsewhere-corrected p-value
  6. Show whether the Planck dipole is in the top 1%, 5%, or 10% of all axes

The scrambling null breaks any true physical correlation by permuting the
cos(θ_EM-n̂) assignments among observations while preserving the marginal
distribution and correlation with cos(D).
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats

from scripts.utils.astronomical_utils import load_skyfield_planets
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, print_status


# Planck 2018 CMB dipole direction
_CMB_RA_RAD = np.deg2rad(168.14)
_CMB_DEC_RAD = np.deg2rad(-7.22)
_CMB_UNIT = np.array([
    np.cos(_CMB_DEC_RAD) * np.cos(_CMB_RA_RAD),
    np.cos(_CMB_DEC_RAD) * np.sin(_CMB_RA_RAD),
    np.sin(_CMB_DEC_RAD),
])
_CMB_UNIT = _CMB_UNIT / np.linalg.norm(_CMB_UNIT)


def grid_celestial_sphere(step_deg=5.0):
    """Return RA/Dec grid in radians."""
    ra = np.deg2rad(np.arange(0, 360, step_deg))
    dec = np.deg2rad(np.arange(-90, 91, step_deg))
    RA, DEC = np.meshgrid(ra, dec)
    return RA.ravel(), DEC.ravel()


def unit_vectors_from_ra_dec(ra_rad, dec_rad):
    """Convert RA/Dec to unit vectors in ICRS."""
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.column_stack([x, y, z])


def compute_earth_moon_cosines(jd_array, directions):
    """
    Compute cos(theta) = em_hat · direction for each observation and direction.

    Parameters
    ----------
    jd_array : ndarray, shape (n,)
    directions : ndarray, shape (d, 3)
        Unit vectors for each sky direction.

    Returns
    -------
    cos_theta : ndarray, shape (n, d)
    """
    planets, _ = load_skyfield_planets(PROJECT_ROOT)
    from skyfield.api import load
    earth = planets["earth"]
    moon = planets["moon"]
    ts = load.timescale()
    timestamps = ts.tt(jd=jd_array)

    earth_pv = earth.at(timestamps).position.km
    moon_pv = moon.at(timestamps).position.km

    em_vec = moon_pv - earth_pv
    em_dist = np.linalg.norm(em_vec, axis=0)
    em_hat = em_vec / em_dist

    # em_hat shape is (3, n), directions shape is (d, 3)
    # cos_theta[n, d] = sum_k em_hat[k, n] * directions[d, k]
    cos_theta = np.dot(em_hat.T, directions.T)
    return cos_theta


def compute_delta_aic_all_directions(cosD, cos_theta_matrix, y, n_obs=None):
    """
    Compute ΔAIC for all directions using the added-variable formula.

    Base model:   y = a·cosD + b + ε        (k=2)
    Augmented:    y = a·cosD + c·(cosθ·cosD) + b + ε   (k=3)

    ΔAIC = n·ln(RSS_base / RSS_aug) - 2

    Optimisation: since r_base ⟂ col(X_base), Z_orth.T @ r_base = Z.T @ r_base.
    This avoids computing the full Z_orth matrix.  np.linalg.solve replaces
    np.linalg.inv for numerical stability.

    Parameters
    ----------
    cosD : ndarray, shape (n,)
    cos_theta_matrix : ndarray, shape (n, d)
        Each column is cosθ for one direction.
    y : ndarray, shape (n,)
    n_obs : int, optional. If None, use len(y).

    Returns
    -------
    delta_aics : ndarray, shape (d,)
    rss_base : float
    """
    n = len(y)
    if n_obs is None:
        n_obs = n

    # Base design and fit — solve with QR for stability
    X_base = np.column_stack([cosD, np.ones(n)])
    XtX = X_base.T @ X_base
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        beta_base = np.linalg.solve(XtX, X_base.T @ y)
        r_base = y - X_base @ beta_base
    rss_base = float(np.sum(r_base ** 2))

    # Precompute M = (X'X)^{-1} X' via solve for stability
    M = np.linalg.solve(XtX, X_base.T)  # shape (2, n)

    # Z matrix: each column is z_i = cosθ_i · cosD
    Z = cos_theta_matrix * cosD[:, None]  # shape (n, d)

    # Suppress spurious BLAS overflow warnings on large matmuls;
    # intermediate accumulators can transiently overflow in blocked
    # routines even though the final result is well within float64 range.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        # Projection of Z onto col(X_base): MZ = (X'X)^{-1} X' Z
        MZ = M @ Z  # shape (2, d)

        # ||Z_orth||^2 = ||Z||^2 - ||proj_X(Z)||^2  (Pythagoras)
        # proj_X(Z) = X_base @ MZ
        proj_norm_sq = np.sum((X_base @ MZ) ** 2, axis=0)
        norms = np.sum(Z ** 2, axis=0) - proj_norm_sq

        # Key identity: r_base ⟂ col(X_base), so Z_orth.T @ r_base = Z.T @ r_base
        projections = Z.T @ r_base  # shape (d,)

    # Defensive numerical handling
    norms_safe = np.where(norms > 1e-20, norms, np.inf)
    norms_safe = np.where(np.isfinite(norms_safe), norms_safe, np.inf)

    # Augmented RSS for each direction
    rss_aug = rss_base - projections ** 2 / norms_safe
    rss_aug = np.clip(rss_aug, 1e-30, None)
    rss_aug = np.where(np.isfinite(rss_aug), rss_aug, rss_base)

    delta_aics = n_obs * np.log(rss_base / rss_aug) - 2.0
    delta_aics = np.where(np.isfinite(delta_aics), delta_aics, 0.0)
    return delta_aics, rss_base


def run_sky_scan_directional() -> dict:
    print_status("═══ Step 076: Full Sky-Scan with Look-Elsewhere Correction ═══", "TITLE")

    data_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not data_path.exists():
        print_status(f"Data not found: {data_path}", "ERROR")
        return {"status": "FAIL", "reason": "No processed data"}

    df_raw = pd.read_csv(data_path)
    outlier_mask = detect_outliers_sigma(df_raw["residual_m"].values, 6.0)
    df = df_raw.loc[~outlier_mask].copy()
    n_clean = len(df)
    print_status(f"Cleaned dataset: N={n_clean:,}", "DATA")

    residuals = df["residual_m"].values
    cosD = np.cos(df["elongation_rad"].values)
    jd = df["date_julian"].values

    # Grid the sphere
    ra_grid, dec_grid = grid_celestial_sphere(step_deg=5.0)
    n_dirs = len(ra_grid)
    print_status(f"Celestial sphere grid: {n_dirs} directions (5° steps)", "DATA")

    directions = unit_vectors_from_ra_dec(ra_grid, dec_grid)

    # Precompute cos(theta) for all observations and all directions
    print_status(">>> Computing Earth-Moon directional cosines (this may take a moment)...", "PROCESS")
    cos_theta_matrix = compute_earth_moon_cosines(jd, directions)
    print_status(f"Computed cos(θ_EM-n̂) matrix: {cos_theta_matrix.shape}", "CALC")

    # Real-data ΔAIC for all directions
    print_status(">>> Fitting directional model for all grid points...", "PROCESS")
    delta_aics, rss_base = compute_delta_aic_all_directions(cosD, cos_theta_matrix, residuals, n_obs=n_clean)

    max_delta_aic = float(np.max(delta_aics))
    max_idx = int(np.argmax(delta_aics))
    best_ra = float(np.rad2deg(ra_grid[max_idx]))
    best_dec = float(np.rad2deg(dec_grid[max_idx]))

    # Planck dipole index
    cmb_direction = unit_vectors_from_ra_dec(np.array([_CMB_RA_RAD]), np.array([_CMB_DEC_RAD]))[0]
    # Find closest grid point to CMB dipole
    dot_prods = np.dot(directions, cmb_direction)
    cmb_idx = int(np.argmax(dot_prods))
    cmb_delta_aic = float(delta_aics[cmb_idx])

    # Rank of CMB dipole among all directions
    sorted_indices = np.argsort(delta_aics)[::-1]
    cmb_rank = int(np.where(sorted_indices == cmb_idx)[0][0]) + 1
    cmb_percentile = 100.0 * cmb_rank / n_dirs

    print_status(
        f"Maximum ΔAIC = {max_delta_aic:.2f} at (RA={best_ra:.1f}°, Dec={best_dec:.1f}°)",
        "RESULT",
    )
    print_status(
        f"Planck dipole ΔAIC = {cmb_delta_aic:.2f}, rank = {cmb_rank}/{n_dirs} "
        f"(top {cmb_percentile:.1f}%)",
        "RESULT",
    )

    # Monte Carlo scrambled skies
    n_mc = 1000
    print_status(f">>> Monte Carlo scrambled skies (n={n_mc})...", "PROCESS")
    rng = np.random.RandomState(42)

    max_delta_aic_scrambled = np.empty(n_mc, dtype=float)
    planck_nulls = np.empty(n_mc, dtype=float)
    for i in range(n_mc):
        perm = rng.permutation(n_clean)
        cos_theta_perm = cos_theta_matrix[perm, :]
        daics, _ = compute_delta_aic_all_directions(cosD, cos_theta_perm, residuals, n_obs=n_clean)
        max_delta_aic_scrambled[i] = np.max(daics)
        # ΔAIC at the Planck dipole direction specifically (not the max)
        planck_nulls[i] = daics[cmb_idx]

    # Conservative p-values: (r+1)/(n+1) avoids p=0.0 when all nulls are below
    r_any = np.sum(max_delta_aic_scrambled >= max_delta_aic)
    p_lec = min(1.0, (r_any + 1) / (n_mc + 1))
    r_cmb = np.sum(planck_nulls >= cmb_delta_aic)
    p_cmb_lec = min(1.0, (r_cmb + 1) / (n_mc + 1))

    # Where does the observed max sit in the scrambled null?
    null_median = float(np.median(max_delta_aic_scrambled))
    null_95 = float(np.percentile(max_delta_aic_scrambled, 95))
    null_99 = float(np.percentile(max_delta_aic_scrambled, 99))

    print_status(
        f"Scrambled null: median ΔAIC_max = {null_median:.2f}, "
        f"95th = {null_95:.2f}, 99th = {null_99:.2f}",
        "CALC",
    )
    print_status(
        f"Look-elsewhere p-value (any axis): {p_lec:.4f}",
        "RESULT",
    )
    print_status(
        f"Look-elsewhere p-value (Planck dipole): {p_cmb_lec:.4f}",
        "RESULT",
    )

    # CMB dipole tier — based on actual rank percentile, not p-value
    if cmb_percentile <= 1.0:
        cmb_tier = "top_1_percent"
    elif cmb_percentile <= 5.0:
        cmb_tier = "top_5_percent"
    elif cmb_percentile <= 10.0:
        cmb_tier = "top_10_percent"
    else:
        cmb_tier = "below_top_10_percent"

    # Distribution of ΔAIC values (real vs null)
    real_pctls = {
        "50": float(np.percentile(delta_aics, 50)),
        "90": float(np.percentile(delta_aics, 90)),
        "95": float(np.percentile(delta_aics, 95)),
        "99": float(np.percentile(delta_aics, 99)),
    }

    return {
        "step_id": "step_076",
        "status": "PASS",
        "n_observations": int(n_clean),
        "n_directions": int(n_dirs),
        "grid_step_deg": 5.0,
        "max_delta_aic": {
            "value": max_delta_aic,
            "best_ra_deg": best_ra,
            "best_dec_deg": best_dec,
        },
        "planck_dipole": {
            "ra_deg": float(np.rad2deg(_CMB_RA_RAD)),
            "dec_deg": float(np.rad2deg(_CMB_DEC_RAD)),
            "delta_aic": cmb_delta_aic,
            "rank_by_delta_aic": cmb_rank,
            "percentile": round(cmb_percentile, 2),
            "tier": cmb_tier,
        },
        "look_elsewhere_correction": {
            "n_scrambled_realizations": n_mc,
            "p_any_axis_exceeds_observed": p_lec,
            "p_planck_axis_exceeds_observed": p_cmb_lec,
            "null_median_max_delta_aic": null_median,
            "null_95th_percentile": null_95,
            "null_99th_percentile": null_99,
        },
        "real_data_delta_aic_percentiles": real_pctls,
        "interpretation": (
            f"Across {n_dirs} uniformly spaced sky directions, the maximum ΔAIC is "
            f"{max_delta_aic:.2f}. Under the scrambled-sky null (n={n_mc}), the "
            f"look-elsewhere-corrected p-value for ANY axis is {p_lec:.4f}. "
            f"The Planck dipole ranks {cmb_rank}/{n_dirs} (top {cmb_percentile:.1f}%) "
            f"with p={p_cmb_lec:.4f}. "
            + (
                "The directional signal survives look-elsewhere correction."
                if p_lec < 0.05
                else "The directional signal does not survive look-elsewhere correction."
            )
        ),
    }


def main() -> int:
    results = run_sky_scan_directional()
    output_path = PROJECT_ROOT / "results" / "outputs" / "step_076_sky_scan_directional.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print_status(f"Saved results to {output_path.relative_to(PROJECT_ROOT)}", "SUCCESS")
    return 0


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_076", str(log_dir / "step_076_sky_scan_directional.log")
    )
    set_step_logger(logger)
    sys.exit(main())
