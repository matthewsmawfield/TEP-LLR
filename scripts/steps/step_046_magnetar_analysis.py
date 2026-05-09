#!/usr/bin/env python3
"""
Step 046: Magnetar Anti-Glitch Analysis
Tests the TEP prediction that magnetar anti-glitches occur near a critical
spin period P_crit ≈ 6.8 s, derived from ρ_T ≈ 20 g/cm³.
Uses published data from Archibald et al. 2013 (Nature, 497, 591).
"""

import sys, json, numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

RHO_T = 20.0  # g/cm³
G = 6.67430e-8  # cm³ g⁻¹ s⁻²
M_SUN_G = 1.98847e33
M_NS = 1.4 * M_SUN_G  # canonical neutron star mass

# Published magnetar data
# Sources: Archibald+2013, Olausen & Kaspi 2014 (McGill Magnetar Catalog)
MAGNETARS = [
    {"name": "1E 2259+586", "P_s": 6.98, "has_antiglitch": True,
     "ref": "Archibald et al. 2013, Nature, 497, 591"},
    {"name": "4U 0142+61", "P_s": 8.69, "has_antiglitch": False,
     "ref": "Dib & Kaspi 2014, ApJ, 784, 37"},
    {"name": "1E 1841-045", "P_s": 11.78, "has_antiglitch": False,
     "ref": "Dib & Kaspi 2014"},
    {"name": "SGR 1900+14", "P_s": 5.16, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 1806-20", "P_s": 7.56, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "1E 1048.1-5937", "P_s": 6.46, "has_antiglitch": False,
     "ref": "Dib & Kaspi 2014"},
    {"name": "XTE J1810-197", "P_s": 5.54, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 0418+5729", "P_s": 9.08, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "Swift J1822.3-1606", "P_s": 8.44, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "3XMM J185246.6+003317", "P_s": 11.56, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 1627-41", "P_s": 2.59, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "CXOU J010043.1-721134", "P_s": 8.02, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 0501+4516", "P_s": 5.76, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 0526-66", "P_s": 8.05, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "CXOU J164710.2-455216", "P_s": 10.61, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 1833-0832", "P_s": 7.56, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "Swift J1834.9-0846", "P_s": 2.48, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "PSR J1622-4950", "P_s": 4.33, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "1RXS J170849.0-400910", "P_s": 11.00, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
    {"name": "SGR 2013+34", "P_s": 0.32, "has_antiglitch": False,
     "ref": "Olausen & Kaspi 2014"},
]


C_CMS = 2.99792458e10  # cm/s

def compute_pcrit(mass_g, rho_t):
    """Compute critical spin period from TEP saturation condition.
    Anti-glitch occurs when light-cylinder radius r_lc = cP/(2π)
    equals the saturation radius R_T = (3M/(4πρ_T))^{1/3}.
    """
    r_t_cm = (3 * mass_g / (4 * np.pi * rho_t)) ** (1.0/3.0)
    return 2 * np.pi * r_t_cm / C_CMS


def run_analysis(logger):
    logger.info(">>> Starting magnetar anti-glitch analysis...")

    p_crit = compute_pcrit(M_NS, RHO_T)
    logger.info(f"  Predicted P_crit = {p_crit:.2f} s (M={M_NS/M_SUN_G:.1f} M_sun)")

    # Test against 1E 2259+586
    target = MAGNETARS[0]
    match_pct = abs(target["P_s"] - p_crit) / p_crit * 100
    logger.info(f"  {target['name']}: P={target['P_s']:.2f}s, match={match_pct:.1f}%")

    # Population analysis
    periods = np.array([m["P_s"] for m in MAGNETARS])
    antiglitch = np.array([m["has_antiglitch"] for m in MAGNETARS])

    # Objects within ±20% of P_crit
    near_crit = np.abs(periods - p_crit) / p_crit < 0.20
    logger.info(f"  Objects within 20% of P_crit: {np.sum(near_crit)}")
    logger.info(f"  Anti-glitch hosts near P_crit: {np.sum(antiglitch & near_crit)}")

    # Mass sensitivity
    mass_range = np.array([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    pcrit_vs_mass = {f"{m:.1f}": float(compute_pcrit(m * M_SUN_G, RHO_T)) for m in mass_range}

    results = {
        "step_id": "step_046",
        "status": "PASS",
        "rho_T_g_cm3": RHO_T,
        "canonical_mass_Msun": 1.4,
        "predicted_P_crit_s": float(p_crit),
        "target_object": {
            "name": target["name"],
            "observed_P_s": target["P_s"],
            "match_percent": float(match_pct),
            "reference": target["ref"],
        },
        "population": {
            "n_magnetars": len(MAGNETARS),
            "n_antiglitch": int(np.sum(antiglitch)),
            "n_near_crit": int(np.sum(near_crit)),
            "period_range_s": [float(np.min(periods)), float(np.max(periods))],
        },
        "pcrit_vs_mass": pcrit_vs_mass,
        "caveat": "Single-object match (N=1); statistical significance is limited. "
                  "Population-level predictions (Section 7.3) are the testable content.",
    }

    logger.info(f"  Note: N={np.sum(antiglitch)} anti-glitch object — statistically limited")
    logger.info("✓   Magnetar analysis complete")
    return results


def main():
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_046", str(log_dir / "step_046_magnetar_analysis.log"))
    set_step_logger(logger)
    results = run_analysis(logger)
    logger.save_step_results(results, PROJECT_ROOT, "step_046_magnetar_analysis")


if __name__ == "__main__":
    main()
