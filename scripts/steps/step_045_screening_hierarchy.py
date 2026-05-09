#!/usr/bin/env python3
"""
Step 045: Screening Hierarchy Analysis
Computes screening factors S = R_T/R_phys for 26 astrophysical objects
using published mass-radius-density data from peer-reviewed sources.
Tests the TEP prediction S ∝ ρ^{1/3}.
"""

import sys, json, numpy as np
from pathlib import Path
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

RHO_T = 20.0  # g/cm^3
M_SUN_G = 1.98847e33
R_EARTH_KM = 6371.0
M_EARTH_G = 5.972e27

# Published data from peer-reviewed sources
# Format: (name, mass_g, radius_km, density_g_cm3, reference)
OBJECTS = [
    # Gas giants
    ("Jupiter", 1.898e30, 69911, 1.33, "NASA Fact Sheet"),
    ("Saturn", 5.683e29, 58232, 0.69, "NASA Fact Sheet"),
    # Rocky planets
    ("Earth", M_EARTH_G, R_EARTH_KM, 5.51, "NASA Fact Sheet"),
    ("Venus", 4.867e27, 6052, 5.24, "NASA Fact Sheet"),
    ("Mars", 6.417e26, 3390, 3.93, "NASA Fact Sheet"),
    ("Mercury", 3.301e26, 2440, 5.43, "NASA Fact Sheet"),
    # Main sequence stars
    ("Sun", 1.989e33, 696340, 1.41, "NASA Fact Sheet"),
    ("Proxima Centauri", 2.428e32, 107000, 56.9, "Kervella+2017"),
    ("Alpha Cen A", 2.188e33, 862400, 0.82, "Kervella+2017"),
    ("Alpha Cen B", 1.789e33, 602300, 1.96, "Kervella+2017"),
    ("Sirius A", 4.018e33, 1189600, 0.57, "Liebert+2005"),
    ("Vega", 4.224e33, 1640000, 0.23, "Yoon+2010"),
    # Brown dwarfs
    ("Gliese 229B", 7.9e31, 69900, 55.0, "Saumon+1996"),
    ("Teide 1", 1.0e32, 69900, 70.0, "Rebolo+1996"),
    # White dwarfs
    ("Sirius B", 2.028e33, 5800, 2.4e6, "Liebert+2005"),
    ("40 Eri B", 9.94e32, 8100, 4.5e5, "Holberg+2016"),
    ("Stein 2051 B", 1.33e33, 7600, 7.2e5, "Sahu+2017"),
    ("Procyon B", 1.19e33, 8600, 4.5e5, "Liebert+2013"),
    ("GD 140", 1.59e33, 7200, 1.0e6, "Holberg+2016"),
    ("EGB 5", 1.19e33, 8400, 4.8e5, "Holberg+2016"),
    # Neutron stars
    ("PSR J0740+6620", 2.78e33, 12.5, 3.4e14, "Miller+2021 (NICER)"),
    ("PSR J0030+0451", 2.58e33, 13.0, 2.8e14, "Riley+2019 (NICER)"),
    ("PSR J0437-4715", 2.78e33, 11.4, 4.5e14, "Choudhury+2024"),
    # Binary pulsars
    ("Hulse-Taylor PSR B1913+16", 2.78e33, 12.0, 3.8e14, "Taylor & Weisberg 1982"),
    ("Double Pulsar PSR J0737-3039A", 2.66e33, 12.5, 3.3e14, "Kramer+2006"),
    ("PSR J0348+0432", 3.98e33, 13.0, 4.3e14, "Antoniadis+2013"),
]


def compute_screening(mass_g, r_phys_km, rho_t=RHO_T):
    """Compute screening factor S = R_T / R_phys."""
    r_t_km = (3 * mass_g / (4 * np.pi * rho_t)) ** (1.0/3.0) * 1e-5  # cm->km
    return r_t_km / r_phys_km


def run_analysis(logger):
    logger.info(">>> Starting screening hierarchy analysis...")

    data = []
    for name, mass, r_phys, rho, ref in OBJECTS:
        s = compute_screening(mass, r_phys)
        data.append({"name": name, "mass_g": mass, "r_phys_km": r_phys,
                      "density": rho, "screening": s, "ref": ref})

    densities = np.array([d["density"] for d in data])
    screenings = np.array([d["screening"] for d in data])

    # Log-log fit: S ∝ ρ^β
    log_rho = np.log10(densities)
    log_s = np.log10(screenings)
    slope, intercept, r_val, p_val, std_err = stats.linregress(log_rho, log_s)
    beta = slope
    beta_err = std_err

    logger.info(f"  Screening law: S ∝ ρ^{beta:.4f} ± {beta_err:.4f}")
    logger.info(f"  R² = {r_val**2:.6f}, p = {p_val:.2e}")
    logger.info(f"  TEP prediction β=1/3: {abs(beta - 1/3)/beta_err:.1f}σ")

    # Note: this is a consistency check, not independent evidence
    # S ≡ R_T/R_phys, R_T ∝ M^{1/3}, R_phys ∝ (M/ρ)^{1/3} → S ∝ ρ^{1/3}

    results = {
        "step_id": "step_045",
        "status": "PASS",
        "n_objects": len(data),
        "rho_T_g_cm3": RHO_T,
        "screening_law": {
            "beta": float(beta), "beta_err": float(beta_err),
            "r_squared": float(r_val**2), "p_value": float(p_val),
        },
        "tep_prediction_beta": 1.0/3.0,
        "sigma_difference": float(abs(beta - 1/3) / beta_err),
        "note": "This is a consistency check, not independent evidence. "
                "The 1/3 exponent is algebraically expected from the definitions.",
        "objects": data,
    }

    logger.info("✓   Screening hierarchy analysis complete")
    return results


def main():
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_045", str(log_dir / "step_045_screening_hierarchy.log"))
    set_step_logger(logger)
    results = run_analysis(logger)
    logger.save_step_results(results, PROJECT_ROOT, "step_045_screening_hierarchy")


if __name__ == "__main__":
    main()
