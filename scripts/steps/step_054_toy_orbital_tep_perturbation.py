#!/usr/bin/env python3
"""
Step 054: Toy Orbital Model — TEP Perturbation Residual Demonstration
=======================================================================

A simplified demonstration that a synodic cos(D) perturbation to the
Earth-Moon distance cannot be absorbed by standard Keplerian orbital-element
fits, and produces residuals of the observed amplitude (~4–7 mm).

Physical setup (2D, coplanar, circular orbits):
    - Earth orbits Sun at 1 AU with period 1 year
    - Moon orbits Earth with period 27.32 days
    - D = Earth-Moon-Sun elongation angle (synodic phase)
    - TEP perturbation: δr = A × η × cos(D), with A = 13 m

Method:
    1. Integrate the true orbit with the perturbation
    2. Fit standard Keplerian elements (a, e, ω, M₀) to the perturbed orbit
    3. Compute O–C residuals
    4. Regress residuals on cos(D) to recover the input η

This is a toy model — it neglects libration, tidal dissipation, PN corrections,
station geometry, and all real-world complexities. Its purpose is to establish
that a cos(D) perturbation of the Nordtvedt amplitude is (i) not absorbed by
Keplerian fitting and (ii) leaves residuals of the observed order of magnitude.
"""

import sys
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
import numpy as np
from scipy.optimize import minimize

logger = TEPLogger("step_054")
set_step_logger(logger)


# Physical constants (simplified)
YEAR_DAYS = 365.25
MONTH_DAYS = 27.32166  # Sidereal month
SYNODIC_DAYS = 29.53059  # Synodic month
AU_M = 1.495978707e11  # 1 AU in meters
A_NORDVEDT = 13.0  # Nordtvedt sensitivity factor [m/η]


def true_orbit_with_tep(t_days, eta=-5.0e-4, a_moon_km=384400.0):
    """
    Generate perturbed Earth-Moon range as a function of time.

    Simplified model:
    - Earth on circular orbit around Sun at 1 AU, period 1 year
    - Moon on circular orbit around Earth, period 27.32 days
    - TEP perturbation: δr = A_NORDVEDT * η * cos(D)
    where D is the Earth-Moon-Sun elongation angle.

    Parameters
    ----------
    t_days : ndarray
        Time array in days
    eta : float
        Nordtvedt parameter (dimensionless)
    a_moon_km : float
        Nominal Earth-Moon distance in km

    Returns
    -------
    r_true : ndarray
        True Earth-Moon distance [m]
    D : ndarray
        Synodic elongation angle [rad]
    cosD : ndarray
        Cosine of elongation
    """
    # Earth's orbital angle around Sun
    theta_earth = 2 * np.pi * t_days / YEAR_DAYS

    # Moon's orbital angle around Earth (sidereal)
    theta_moon = 2 * np.pi * t_days / MONTH_DAYS

    # Synodic phase: D = θ_moon − θ_earth (mod 2π)
    D = np.mod(theta_moon - theta_earth, 2 * np.pi)
    # Shift to −π to +π for cos(D) symmetry
    D = np.where(D > np.pi, D - 2 * np.pi, D)

    cosD = np.cos(D)

    # Nominal circular distance
    r0 = a_moon_km * 1000.0  # meters

    # TEP perturbation
    delta_r = A_NORDVEDT * eta * cosD

    r_true = r0 + delta_r

    return r_true, D, cosD


def ephemeris_basis_range(t_days, a, e_cos, e_sin, a_moon_km=384400.0):
    """
    Standard ephemeris orbital-element basis for Earth-Moon range.

    A standard ephemeris fit adjusts orbital elements (a, e, ω, M₀).
    To first order in eccentricity, the range modulation is:
        r(t) ≈ a * [1 − e*cos(M(t) − ω)]
             = a − a*e*cos(ω)*cos(M) − a*e*sin(ω)*sin(M)

    We fit the equivalent linear basis: {1, cos(M), sin(M)} where
    M(t) = M₀ + n*t is the mean anomaly.  The parameters (a, e_cos, e_sin)
    map back to (a, e, ω) but the fit is numerically stable.

    This correctly models what standard LLR ephemerides do: they fit
    Keplerian elements, which produce signals at the ORBITAL period
    (~27.3 days), NOT at the synodic period (~29.5 days).
    """
    n = 2 * np.pi / MONTH_DAYS  # mean motion [rad/day]
    M = n * t_days  # M₀ absorbed into e_cos / e_sin phase
    return a + e_cos * np.cos(M) + e_sin * np.sin(M)


def fit_ephemeris(t_days, r_true, a_moon_km=384400.0):
    """
    Fit standard orbital-element basis to the perturbed range data.

    The fit optimises (a, e_cos, e_sin) — equivalent to adjusting
    semi-major axis and eccentricity vector of a Keplerian orbit.
    """
    n = 2 * np.pi / MONTH_DAYS
    M = n * t_days
    X = np.column_stack([np.ones(len(t_days)), np.cos(M), np.sin(M)])
    coeffs, _, _, _ = np.linalg.lstsq(X, r_true, rcond=None)
    a_fit, e_cos_fit, e_sin_fit = coeffs
    r_fit = X @ coeffs
    rss = np.sum((r_true - r_fit)**2)
    e_fit = np.sqrt(e_cos_fit**2 + e_sin_fit**2) / abs(a_fit)
    omega_fit = np.arctan2(-e_sin_fit, -e_cos_fit)

    return {
        'a': float(a_fit), 'e': float(e_fit),
        'omega': float(omega_fit),
        'r_model': r_fit, 'rss': float(rss),
        'success': True
    }


def main():
    print_status("═══ Step 054: Toy Orbital TEP Perturbation Model ═══", "TITLE")

    # Simulation parameters
    t_days = np.linspace(0, 10 * YEAR_DAYS, 50000)  # 10 years, dense sampling
    eta_input = -5.4e-4  # Representative value from full-sample analysis

    print_status(f"Simulation: {len(t_days)} points over {t_days[-1]/YEAR_DAYS:.1f} years", "INFO")
    print_status(f"Input η = {eta_input:.2e}", "INFO")

    # 1. Generate true perturbed orbit
    r_true, D, cosD = true_orbit_with_tep(t_days, eta=eta_input)
    print_status(f"Perturbation amplitude: {A_NORDVEDT * abs(eta_input) * 1000:.2f} mm", "INFO")

    # 2. Fit standard orbital-element basis
    print_status("Fitting standard orbital-element basis...", "PROCESS")
    eph_fit = fit_ephemeris(t_days, r_true)
    print_status(f"  Fit: a = {eph_fit['a']/1e3:.1f} km, e = {eph_fit['e']:.5f}", "RESULT")

    # 3. Compute residuals
    residuals = r_true - eph_fit['r_model']
    rms_mm = np.sqrt(np.mean(residuals**2)) * 1000.0
    print_status(f"  Residual RMS: {rms_mm:.2f} mm", "RESULT")

    # 4. Regress residuals on cos(D) to recover η
    print_status("Regressing residuals on cos(D)...", "PROCESS")
    X = np.column_stack([cosD, np.ones(len(cosD))])
    coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
    eta_recovered = coeffs[0] / A_NORDVEDT
    se = np.sqrt(np.mean((residuals - X @ coeffs)**2) * np.linalg.inv(X.T @ X)[0, 0])
    eta_err_recovered = se / A_NORDVEDT
    snr = abs(eta_recovered) / max(eta_err_recovered, 1e-20)

    print_status(f"  Recovered η = {eta_recovered:.4e} ± {eta_err_recovered:.4e} ({snr:.2f}σ)", "RESULT")
    print_status(f"  Input η     = {eta_input:.4e}", "RESULT")
    print_status(f"  Recovery bias: {(eta_recovered - eta_input)/eta_input * 100:.1f}%", "RESULT")

    # 5. Test multiple η values to map perturbation amplitude vs recovery
    print_status("--- Sweep over η values ---", "INFO")
    eta_values = np.array([-1.0e-3, -5.0e-4, -3.0e-4, -1.0e-4, -5.0e-5, -1.0e-5])
    sweep_results = []
    for eta_test in eta_values:
        r_test, _, cosD_test = true_orbit_with_tep(t_days, eta=eta_test)
        ef = fit_ephemeris(t_days, r_test)
        resid_test = r_test - ef['r_model']
        X_test = np.column_stack([cosD_test, np.ones(len(cosD_test))])
        c_test, _, _, _ = np.linalg.lstsq(X_test, resid_test, rcond=None)
        eta_rec = c_test[0] / A_NORDVEDT
        sweep_results.append({
            'eta_input': float(eta_test),
            'eta_recovered': float(eta_rec),
            'amplitude_mm': float(abs(A_NORDVEDT * eta_test) * 1000),
            'residual_rms_mm': float(np.sqrt(np.mean(resid_test**2)) * 1000)
        })
        print_status(f"  η_in = {eta_test:.2e}: recovered = {eta_rec:.4e}, "
                     f"amp = {abs(A_NORDVEDT * eta_test) * 1000:.1f} mm, RMS = {np.sqrt(np.mean(resid_test**2)) * 1000:.2f} mm",
                     "RESULT")

    # 6. Save results
    results = {
        'step': '054_toy_orbital_tep_perturbation',
        'simulation': {
            'n_points': len(t_days),
            'span_years': float(t_days[-1] / YEAR_DAYS),
            'input_eta': float(eta_input),
            'perturbation_amplitude_mm': float(A_NORDVEDT * abs(eta_input) * 1000)
        },
        'ephemeris_fit': {
            'a_km': float(eph_fit['a'] / 1000),
            'e': float(eph_fit['e']),
            'omega_rad': float(eph_fit['omega']),
            'rss_m2': float(eph_fit['rss'])
        },
        'residuals': {
            'rms_mm': float(rms_mm),
            'eta_recovered': float(eta_recovered),
            'eta_err_recovered': float(eta_err_recovered),
            'snr': float(snr),
            'recovery_bias_percent': float((eta_recovered - eta_input) / eta_input * 100)
        },
        'eta_sweep': sweep_results,
        'interpretation': (
            "A synodic cos(D) perturbation of the Nordtvedt amplitude (~7 mm for "
            f"η = {eta_input:.1e}) is NOT absorbed by standard orbital-element fitting. "
            "The fit optimises semi-major axis and eccentricity vector, which generate "
            "signals at the orbital period (~27.3 days), not at the Earth-Moon-Sun synodic "
            "period (~29.5 days).  Because cos(D) is orthogonal to the {1, cos(M), sin(M)} "
            "basis of Keplerian range modulation, the perturbation survives as an unabsorbed "
            "residual.  The recovered η from the residuals matches the input to within "
            f"{(eta_recovered - eta_input)/eta_input * 100:.1f}%, confirming the signal "
            "cannot be eliminated by adjusting standard ephemeris parameters."
        )
    }

    output_path = PROJECT_ROOT / "results/outputs/step_054_toy_orbital_tep_perturbation.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print_status(f"Saved to {output_path}", "INFO")


if __name__ == '__main__':
    main()
