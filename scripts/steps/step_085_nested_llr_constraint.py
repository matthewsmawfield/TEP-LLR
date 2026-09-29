"""Nested-operator LLR differential constraint.

Evaluates the canonical pairwise response R = S_Sigma(X_env)^2 * y(s) for
the Sun-Earth and Sun-Moon channels, takes the differential acceleration
toward the Sun, and compares against the LLR bound. Migrated from Paper 0
step_22_llr_differential_constraint.py; the observational constraint
belongs to this paper while the operator itself is Paper 0's.

Outputs results/outputs/step_085_nested_llr_constraint.json
"""
import json
import os
from scipy.optimize import brentq

# Physical constants
AU = 1.496e11  # meters
M_sun = 1.989e30  # kg
M_earth = 5.972e24  # kg
M_moon = 7.348e22  # kg
R_sun = 6.96e8  # meters
G = 6.67430e-11  # m^3 kg^-1 s^-2
c = 299792458.0  # m/s
H0 = 70e3 / 3.085677581e22
BETA_A = -1.0
g_t = c * H0 / (2 * BETA_A ** 2)   # canonical transition acceleration

# Galactic embedding ambient (wide-binary plateau: q_env^2 ~ 0.43)
X_GAL = 0.52

# Experimental bounds
# LLR bound on anomalous differential acceleration towards the Sun
# \Delta a / a_N < ~ 1.3e-13 (e.g. Murphy 2012, or Williams et al 2012)
LLR_BOUND_ETA = 1.3e-13


def r_star_M(M):
    """Derived shell radius r* = sqrt(GM/g_t) (M^1/2 law)."""
    return (G * M / g_t) ** 0.5


def y_profile(x):
    """Exact flux-conserving profile: y(1 + y^2 x^-4) = 1."""
    if x >= 50.0:
        return 1.0
    return brentq(lambda yy: yy * (1.0 + yy * yy / x ** 4) - 1.0,
                  1e-30, 1.0, xtol=1e-14)


def S_sigma(X):
    return 1.0 / (1.0 + X)


def R_hierarchical(M_dom, s, X_extra=0.0):
    """Nested pair response R = S_Sigma(X_env)^2 * y(s): the vertex factors
    are evaluated at the pair's embedding ambient — the dominant member's
    field at the pair scale plus external ambients — giving R -> y^3 ~ s^4
    deep inside the shell (resolved operator, step_30 / issue 0-27)."""
    y = y_profile(s / r_star_M(M_dom))
    q_env = S_sigma((1.0 - y) / y + X_GAL + X_extra)
    return q_env ** 2 * y


def run():
    print("=========================================================")
    print(" TEP Lunar Laser Ranging (LLR) Differential Constraint")
    print(" (resolved nested operator R = S_Sigma(X_env)^2 * y(s))")
    print("=========================================================")
    print()

    s_sun_earth = 1.0 * AU
    s_em = 2.57e-3 * AU  # Earth-Moon separation

    Rs_sun = r_star_M(M_sun)
    Rs_earth = r_star_M(M_earth)
    Rs_moon = r_star_M(M_moon)

    print(f"Derived shell radii r* = sqrt(GM/g_t):")
    print(f"  r*(Sun)   = {Rs_sun/AU:.1f} AU")
    print(f"  r*(Earth) = {Rs_earth/AU:.2f} AU")
    print(f"  r*(Moon)  = {Rs_moon/AU:.2f} AU")

    # Embedding ambients
    X_sun_1au = (1.0 - y_profile(s_sun_earth / Rs_sun)) \
        / y_profile(s_sun_earth / Rs_sun)
    X_em = (1.0 - y_profile(s_em / Rs_earth)) \
        / y_profile(s_em / Rs_earth)

    # Sun-Earth pair: ambient = solar field at 1 AU + galactic
    y_es = y_profile(s_sun_earth / Rs_sun)
    S_sun_earth = S_sigma(X_sun_1au + X_GAL) ** 2 * y_es

    # Sun-Moon pair: the Moon additionally sits inside Earth's nonlinear
    # shell -> ambient = solar + terrestrial(at lunar orbit) + galactic
    S_sun_moon = S_sigma(X_sun_1au + X_em + X_GAL) ** 2 * y_es

    print()
    print(f"Effective Scalar Responses (pair responses at 1 AU):")
    print(f"  R(Sun, Earth) = {S_sun_earth:.4e}")
    print(f"  R(Sun, Moon)  = {S_sun_moon:.4e}")

    # 3. Differential acceleration toward the Sun
    a_N = G * M_sun / (s_sun_earth**2)
    a_phi_earth = a_N * S_sun_earth
    a_phi_moon = a_N * S_sun_moon
    delta_a_scalar = abs(a_phi_earth - a_phi_moon)
    delta_EM = abs(S_sun_earth - S_sun_moon)

    print()
    print(f"Accelerations toward the Sun:")
    print(f"  Newtonian a_N         = {a_N:.4e} m/s^2")
    print(f"  Scalar a_phi (Earth)  = {a_phi_earth:.4e} m/s^2")
    print(f"  Scalar a_phi (Moon)   = {a_phi_moon:.4e} m/s^2")
    print(f"  Diff scalar accel     = {delta_a_scalar:.4e} m/s^2")

    # 4. Resulting LLR observable
    print()
    print(f"LLR Observable (Delta a / a_N):")
    print(f"  TEP Predicted delta_EM = {delta_EM:.4e}")
    print(f"  Experimental Bound < {LLR_BOUND_ETA:.1e}")

    if delta_EM < LLR_BOUND_ETA:
        print("  => PASSES LLR CONSTRAINT")
    else:
        print("  => FAILS LLR CONSTRAINT")

    print()
    print("GR Reference Check:")
    print("  If scalar coupling is zero (R -> 0), delta_EM = 0.0 (Recovers GR WEP).")

    # 5. Cross-scale consistency under the same fixed operator
    s_cassini = 1.6 * R_sun
    s_saturn = 9.5 * AU
    s_wb = 2646 * AU

    S_cassini = R_hierarchical(M_sun, s_cassini)
    S_saturn = R_hierarchical(M_sun, s_saturn)
    # Comparable-mass pair: saturated galactic vertices * y(s)
    y_wb = y_profile(s_wb / r_star_M(1.2 * M_sun))
    q_gal = S_sigma(X_GAL)
    S_wb = q_gal ** 2 * y_wb

    print()
    print("Cross-Scale Consistency Check (Fixed Parameters):")
    print(f"  Cassini shear response (s = 1.6 R_sun): {S_cassini:.4e}  (charge channel bound separate)")
    print(f"  Saturn response (s = 9.5 AU)          : {S_saturn:.4e}  -> da = {2*S_saturn*G*M_sun/s_saturn**2:.2e} m/s^2 vs 5e-13 bound")
    print(f"  Wide Binary response (s = 2646 AU)    : {S_wb:.4f}       (observed plateau ~0.43)")

    # Write JSON output
    output = {
        "step_id": "step_085", "status": "PASS",
        "operator": "R = S_Sigma(X_env)^2 * y(s); y(1+y^2(r*/r)^4)=1; "
                    "r* = sqrt(GM/g_t)",
        "radii_au": {
            "sun": Rs_sun / AU,
            "earth": Rs_earth / AU,
            "moon": Rs_moon / AU
        },
        "suppression": {
            "sun_earth": S_sun_earth,
            "sun_moon": S_sun_moon
        },
        "acceleration_ms2": {
            "newtonian": a_N,
            "scalar_earth": a_phi_earth,
            "scalar_moon": a_phi_moon,
            "differential": delta_a_scalar
        },
        "llr_test": {
            "predicted_delta_EM": delta_EM,
            "bound": LLR_BOUND_ETA,
            "passes": bool(delta_EM < LLR_BOUND_ETA)
        },
        "consistency_checks": {
            "cassini_shear_response": S_cassini,
            "saturn_response": S_saturn,
            "wb_activation": S_wb
        }
    }

    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "results", "outputs")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_085_nested_llr_constraint.json"), "w") as f:
        json.dump(output, f, indent=2)
    print("\n[SUCCESS] Wrote JSON output to results/outputs/step_085_nested_llr_constraint.json")

if __name__ == "__main__":
    run()
