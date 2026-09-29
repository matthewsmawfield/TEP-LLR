#!/usr/bin/env python3
"""Step 084: finite-body matched-profile Earth–Moon–Sun bridge diagnostic.

Each uniform-density body's radial kinetic flux relation is solved separately;
its gradient is then combined with the other body's gradient and a locally
uniform solar background. This matched-profile ansatz does not solve the
nonlinear three-source PDE, compute a mutual force, or fit an LLR range model.
The baseline and two-branch kinetic coefficients reproduce Paper 0 step_56's
specified constitutive laws; their constants are pinned locally for standalone
Paper 17 reproducibility.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.utils.llr_constants import C_M_S, GM_EARTH, GM_MOON, GM_SUN, MOON_RADIUS_M

STEP_NUM = "084"
AU_M = 1.495978707e11
SEPARATION_M = 3.84399e8
EARTH_RADIUS_M = 6371e3
SUN_RADIUS_M = 6.96e8
H0_S_INV = 70e3 / 3.085677581e22
G_T = C_M_S * H0_S_INV / 2
K = 16.03
EPS_REG = 1e-8


def kinetic_coefficient(q, sector):
    if sector == "baseline":
        return 1.0 + q * q
    if sector == "two_branch":
        return EPS_REG + K * q / math.sqrt(2) + q * q
    raise ValueError(f"unknown kinetic sector: {sector}")


def radial_gradient(gm, radius, distance, sector):
    if distance < 0:
        raise ValueError("radial distance must be non-negative")
    if distance == 0:
        return 0.0
    enclosed_gm = gm * min(1.0, (distance / radius) ** 3)
    target = enclosed_gm / (G_T * distance ** 2)
    return brentq(lambda q: kinetic_coefficient(q, sector) * q - target,
                  0.0, max(1.0, target), xtol=1e-12)


def bridge(sector, angle_deg):
    solar = radial_gradient(GM_SUN, SUN_RADIUS_M, AU_M, sector)
    angle = math.radians(angle_deg)

    def axial(x):
        earth = radial_gradient(GM_EARTH, EARTH_RADIUS_M, x, sector)
        moon = radial_gradient(GM_MOON, MOON_RADIUS_M, SEPARATION_M - x, sector)
        return -earth + moon + solar * math.cos(angle)

    def field(x):
        return math.hypot(axial(x), solar * math.sin(angle))

    lo, hi = EARTH_RADIUS_M, SEPARATION_M - MOON_RADIUS_M
    if abs(math.sin(angle)) < 1e-12 and axial(lo) * axial(hi) < 0:
        x_min = brentq(axial, lo, hi, xtol=1e-4)
        q_min = field(x_min)
    else:
        grid = np.linspace(lo, hi, 401)
        values = np.array([field(x) for x in grid])
        candidates = [(float(values[0]), float(lo)), (float(values[-1]), float(hi))]
        for i in range(1, len(grid) - 1):
            if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
                fit = minimize_scalar(field, bounds=(grid[i - 1], grid[i + 1]),
                                      method="bounded", options={"xatol": 1.0})
                candidates.append((float(fit.fun), float(fit.x)))
        q_min, x_min = min(candidates)
    return {"solar_gradient_over_gt": solar,
            "bridge_min_gradient_over_gt": q_min,
            "bridge_min_distance_from_earth_km": x_min / 1000,
            "bridge_min_P_X": kinetic_coefficient(q_min, sector),
            "solar_P_X": kinetic_coefficient(solar, sector)}


def main():
    if not 0 < EARTH_RADIUS_M < SEPARATION_M - MOON_RADIUS_M:
        raise ValueError("body surfaces overlap")
    rows = {sector: {str(angle): bridge(sector, angle)
                     for angle in (0, 45, 90, 135, 180)}
            for sector in ("baseline", "two_branch")}
    q_critical = brentq(lambda q: kinetic_coefficient(q, "two_branch") - 1,
                        0, 1)
    angular_half_width = math.degrees(math.asin(
        q_critical / rows["two_branch"]["0"]["solar_gradient_over_gt"]))
    fit_path = ROOT / "results/outputs/step_056_dynamical_integrator_eta_refit.json"
    fit = json.loads(fit_path.read_text())
    eta = fit["inpop19a"]["linearized_integrator_eta"]
    rstar_earth = math.sqrt(GM_EARTH / G_T)
    uniform_cells = math.ceil(8 * SEPARATION_M / MOON_RADIUS_M)
    result = {
        "step_id": f"step_{STEP_NUM}", "status": "PASS",
        "scope": "finite uniform-density radial sources, matched-gradient ansatz; not a nonlinear three-body solve or an LLR range fit",
        "inputs": {"g_t_m_s2": G_T, "GM_earth_m3_s2": GM_EARTH,
                   "GM_moon_m3_s2": GM_MOON, "GM_sun_m3_s2": GM_SUN,
                   "earth_radius_m": EARTH_RADIUS_M, "moon_radius_m": MOON_RADIUS_M,
                   "earth_moon_separation_m": SEPARATION_M,
                   "solar_distance_m": AU_M, "kinetic_two_branch_k": K,
                   "kinetic_regularizer": EPS_REG},
        "resolution_comparison": {
            "rstar_earth_m": rstar_earth,
            "earth_moon_over_rstar": SEPARATION_M / rstar_earth,
            "earth_radius_over_separation": EARTH_RADIUS_M / SEPARATION_M,
            "moon_radius_over_separation": MOON_RADIUS_M / SEPARATION_M,
            "uniform_grid_cells_per_axis_for_four_cells_per_moon_radius_over_2d": uniform_cells,
            "uniform_grid_cells_3d_at_that_spacing": uniform_cells ** 3},
        "orientation_degrees_from_earth_moon_axis": [0, 45, 90, 135, 180],
        "sectors": rows,
        "matched_ansatz_P_X_lt_1_angular_half_width_degrees": angular_half_width,
        "matched_ansatz_P_X_lt_1_critical_gradient_over_gt": q_critical,
        "llr_context_not_a_model_prediction": {
            "inpop19a_eta_resid_fit": eta["eta"],
            "inpop19a_eta_resid_fit_error": eta["eta_error"],
            "inpop19a_eta_resid_source": "results/outputs/step_056_dynamical_integrator_eta_refit.json",
            "fit_scope": fit["scope"]},
        "limitations": (
            "The nonlinear PDE does not permit superposition of spherical "
            "solutions. Solar field is held uniform over the pair; each "
            "body is uniform-density and its charge is assumed proportional "
            "to mass. No scalar self-energy sensitivities, three-dimensional "
            "mutual force, lunar orbit perturbation, ephemeris reduction, "
            "or disformal two-way light time is calculated. Local P_X and "
            "gradient minima cannot be interpreted as LLR range anomalies."
        )}
    path = ROOT / "results/outputs/step_084_earth_moon_solar_bridge.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {path}")
    for sector, scan in rows.items():
        print(sector, [(a, round(r["bridge_min_gradient_over_gt"], 3))
                       for a, r in scan.items()])


if __name__ == "__main__":
    main()
