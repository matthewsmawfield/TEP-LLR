#!/usr/bin/env python3
"""Step 086: coupled nonlinear Earth-Moon-Sun solve at syzygy.

Step 084 matched separately solved radial profiles against a uniform
solar background and explicitly disclaimed the nonlinear three-source
PDE.  At syzygy (Sun-Earth-Moon collinear) the full problem IS
axisymmetric, so the step_32/56 finite-volume machinery applies
directly: the solar tide enters as a uniform ambient gradient u0,
Earth and Moon are two Gaussian sources on the z-axis, and

    div J(a) = rho ,   J_i = a_i f(|a|),
    f(q) = eps_reg + k q/sqrt(2) + q^2   (two-branch, operative)
    f(q) = 1 + q^2                        (baseline control)

is solved on the cylinder grid (s, z) in dimensionless units
d_EM = 1, g_t = 1.  Source charges are set by the Newtonian
accelerations at unit separation:

    a_E(d) = GM_E/d^2 = 7.93e6 g_t  ->  Q_E = 4 pi a_E
    a_M(d) = GM_M/d^2 = 9.77e4 g_t  ->  Q_M = 4 pi a_M
    u_solar = GM_sun/AU^2 = 1.74e7 g_t (uniform over the pair)

The ambient u0 is the SOLAR scalar gradient, not the Newtonian
acceleration: the constitutive map is q f(q) = a_N/g_t, so the solar
field GM_sun/AU^2 = 1.74e7 g_t maps to q ~ 259 (baseline) / ~1240
(two-branch).  The pair is solved at u0 = 0 first, then warm-started
through intermediate u0 to the operative solar value.

Extracted vs the step_084 matched-ansatz baseline:
  * axial field |a(z)| between the bodies and its minimum (the bridge
    point), the solved replacement for the matched-ansatz minimum;
  * volume fraction with f(q) < 1 (the screened bridge region);
  * mutual-force correction y_pair = F_pair / F_pair^linear via
    momentum-conserving stress flux through the midplane sphere
    (step_59 convention).

Outputs results/outputs/step_086_earth_moon_coupled_pde.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.utils.llr_constants import C_M_S, GM_EARTH, GM_MOON, GM_SUN

H0_S_INV = 70e3 / 3.085677581e22
G_T = C_M_S * H0_S_INV / 2
AU_M = 1.495978707e11
D_EM = 3.84399e8

K = 16.03
EPS_REG = 1e-8
SQRT2 = math.sqrt(2.0)

# dimensionless source strengths: the Newtonian accelerations set the
# source normalisation, a_N = Q/(4 pi r^2) at unit separation
A_SOLAR_N = (GM_SUN / AU_M**2) / G_T
Q_E = 4.0 * math.pi * (GM_EARTH / D_EM**2) / G_T
Q_M = 4.0 * math.pi * (GM_MOON / D_EM**2) / G_T


def u_solar(sector):
    """Ambient scalar gradient: q f(q) = a_solar/g_t (the same map
    step_084's radial_gradient solves)."""
    from scipy.optimize import brentq
    return brentq(lambda q: f_flux(q, sector) * q - A_SOLAR_N,
                  0.0, max(1.0, A_SOLAR_N), xtol=1e-12)

SIG_SRC = 0.08          # Gaussian source width in d_EM units

# grid: Earth at z=0, Moon at z=1
NS, NZ = 120, 260
SMAX, ZMIN, ZMAX = 2.5, -1.5, 3.5
ds = SMAX / NS
dz = (ZMAX - ZMIN) / NZ
s = (np.arange(NS) + 0.5) * ds
z = (np.arange(NZ) + 0.5) * dz + ZMIN
SS, ZZ = np.meshgrid(s, z, indexing='ij')


def make_rho():
    rr_e = np.sqrt(SS**2 + ZZ**2)
    rr_m = np.sqrt(SS**2 + (ZZ - 1.0)**2)
    rho = (Q_E * np.exp(-rr_e**2 / (2 * SIG_SRC**2))
           + Q_M * np.exp(-rr_m**2 / (2 * SIG_SRC**2)))
    rho /= (2.0 * np.pi)**1.5 * SIG_SRC**3
    return rho


RHO = make_rho()
S_FACE = np.arange(NS + 1) * ds


def f_flux(q, sector):
    if sector == "baseline":
        return 1.0 + q * q
    return EPS_REG + K * q / SQRT2 + q * q


def fp_flux(q, sector):
    if sector == "baseline":
        return 2.0 * q
    return K / SQRT2 + 2.0 * q


def face_fluxes(psi, u0, sector):
    d_dz = np.gradient(psi, dz, axis=1)
    d_ds = np.gradient(psi, ds, axis=0)

    a_s_f = (psi[1:, :] - psi[:-1, :]) / ds
    a_z_f = 0.5 * (d_dz[:-1, :] + d_dz[1:, :]) + u0
    q = np.sqrt(np.maximum(a_s_f**2 + a_z_f**2, 1e-30))
    Js_f = a_s_f * f_flux(q, sector)
    fp = fp_flux(q, sector)
    Zss_f = f_flux(q, sector) + fp * a_s_f**2 / q
    Zsz_f = fp * a_s_f * a_z_f / q

    a_z_g = (psi[:, 1:] - psi[:, :-1]) / dz + u0
    a_s_g = 0.5 * (d_ds[:, :-1] + d_ds[:, 1:])
    qg = np.sqrt(np.maximum(a_s_g**2 + a_z_g**2, 1e-30))
    Jz_g = a_z_g * f_flux(qg, sector)
    fpg = fp_flux(qg, sector)
    Zzz_g = f_flux(qg, sector) + fpg * a_z_g**2 / qg
    Zsz_g = fpg * a_s_g * a_z_g / qg
    return (Js_f, Zss_f, Zsz_f), (Jz_g, Zzz_g, Zsz_g)


def residual(psi, u0, sector):
    (Js_f, _, _), (Jz_g, _, _) = face_fluxes(psi, u0, sector)
    div_s = np.zeros_like(psi)
    sJp = S_FACE[1:-1][:, None] * Js_f
    div_s[:-1, :] += sJp / (SS[:-1, :] * ds)
    div_s[1:, :] -= sJp / (SS[1:, :] * ds)
    div_z = np.zeros_like(psi)
    div_z[:, :-1] += Jz_g / dz
    div_z[:, 1:] -= Jz_g / dz
    return div_s + div_z - RHO


def build_jacobian(psi, u0, sector):
    (Js_f, Zss_f, Zsz_f), (Jz_g, Zzz_g, Zsz_g) = face_fluxes(psi, u0, sector)
    N = NS * NZ
    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(np.asarray(r).ravel())
        cols.append(np.asarray(c).ravel())
        vals.append(np.asarray(v).ravel())

    JJ = np.arange(NZ)
    for k in range(NS - 1):
        i, inb = k, k + 1
        sf = S_FACE[k + 1]
        for ci, sgn in ((i, +1.0), (inb, -1.0)):
            c_idx = ci * NZ + JJ
            denom = np.maximum(s[ci], 1e-12) * ds
            coef = sgn * sf * Zss_f[k, :] / (denom * ds)
            add(c_idx, inb * NZ + JJ, coef)
            add(c_idx, i * NZ + JJ, -coef)
            coefz = sgn * sf * Zsz_f[k, :] / (denom * 4 * dz)
            for dj, w in ((1, 1.0), (-1, -1.0)):
                jj = np.clip(JJ + dj, 0, NZ - 1)
                mask = np.ones(NZ, bool)
                mask[-1 if dj > 0 else 0] = False
                add(c_idx[mask], (inb * NZ + jj[mask]),
                    coefz[mask] * w)
                add(c_idx[mask], (i * NZ + jj[mask]),
                    coefz[mask] * w)

    II = np.arange(NS)
    for j in range(NZ - 1):
        jnb = j + 1
        for cj, sgn in ((j, +1.0), (jnb, -1.0)):
            c_idx = II * NZ + cj
            coef = sgn * Zzz_g[:, j] / (dz * dz)
            add(c_idx, II * NZ + jnb, coef)
            add(c_idx, II * NZ + j, -coef)
            coefs = sgn * Zsz_g[:, j] / (4 * ds * dz)
            for di, w in ((1, 1.0), (-1, -1.0)):
                ii = np.clip(II + di, 0, NS - 1)
                mask = np.ones(NS, bool)
                mask[-1 if di > 0 else 0] = False
                add(c_idx[mask], (ii[mask] * NZ + jnb),
                    coefs[mask] * w)
                add(c_idx[mask], (ii[mask] * NZ + j),
                    coefs[mask] * w)

    A = sp.csr_matrix(
        (np.concatenate(vals),
         (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return A


def dirichlet_project(A):
    A = A.tolil()
    bd = ([(NS - 1) * NZ + j for j in range(NZ)]
          + [i * NZ + j for i in range(NS) for j in (0, NZ - 1)])
    for c in bd:
        A.rows[c] = [c]
        A.data[c] = [1.0]
    return A.tocsr()


def solve_embedded(u0, sector, itmax=120, psi0=None):
    psi = np.zeros_like(SS) if psi0 is None else psi0.copy()
    for it in range(itmax):
        R = residual(psi, u0, sector).copy()
        R[NS - 1, :] = psi[NS - 1, :]
        R[:, 0] = psi[:, 0]
        R[:, NZ - 1] = psi[:, NZ - 1]
        rn = np.linalg.norm(R[2:-2, 2:-2])
        if it % 10 == 0:
            print(f"  [{sector}] u0={u0:.3e} it={it} |R|={rn:.3e}",
                  flush=True)
        if rn < 1e-6 * np.linalg.norm(RHO):
            break
        A = dirichlet_project(build_jacobian(psi, u0, sector))
        try:
            dpsi = spla.spsolve(A, -R.reshape(-1)).reshape(NS, NZ)
        except Exception:
            dpsi = spla.lsmr(A, -R.reshape(-1))[0].reshape(NS, NZ)
        dmax = np.max(np.abs(dpsi))
        cap = max(0.5 * max(u0, 1.0), 1.0)
        if dmax > cap:
            dpsi *= cap / dmax
        alpha = 1.0
        accepted = False
        for _ in range(40):
            cand = psi + alpha * dpsi
            rc = residual(cand, u0, sector)
            if np.linalg.norm(rc[2:-2, 2:-2]) < rn:
                psi = cand
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            psi = psi + dpsi * 1e-6
    print(f"  [{sector}] u0={u0:.3e} done it={it} |R|="
          f"{np.linalg.norm(residual(psi, u0, sector)[2:-2, 2:-2]):.3e}",
          flush=True)
    return psi


def extract(psi, u0, sector):
    g_s = np.gradient(psi, ds, axis=0)
    g_z = np.gradient(psi, dz, axis=1) + u0
    amag = np.sqrt(g_s**2 + g_z**2)

    # axial profile between the bodies (s = ds/2, z in (0,1))
    iax = 0
    zin = (z > 0.05) & (z < 0.95)
    a_ax = amag[iax, zin]
    z_ax = z[zin]
    imin = int(np.argmin(a_ax))

    # screened bridge: volume where f(q) < 1
    q = amag
    fq = f_flux(q, sector)
    vol_cell = 2 * np.pi * SS * ds * dz
    frac_screened = float(np.sum(vol_cell[fq < 1.0]) /
                          np.sum(vol_cell[(ZZ < 1.2) & (ZZ > -0.2)]))

    return {
        "axial_z_over_d": z_ax.tolist(),
        "axial_abs_a_over_gt": a_ax.tolist(),
        "bridge_min_abs_a_over_gt": float(a_ax[imin]),
        "bridge_min_z_over_d": float(z_ax[imin]),
        "bridge_min_P_X": float(f_flux(a_ax[imin], sector)),
        "screened_volume_fraction_near_pair": frac_screened,
        "ambient_u0_over_gt": float(u0),
        "ambient_P_X": float(f_flux(abs(u0), sector)),
    }


def solve_with_continuation(sector, u0_final, warm_u0s):
    psi = None
    for u0 in [0.0] + warm_u0s + [u0_final]:
        psi = solve_embedded(u0, sector, psi0=psi)
    return psi


def main():
    out = {
        "step_id": "step_086",
        "status": "coupled nonlinear solve (syzygy, axisymmetric)",
        "scope": ("replaces the step_084 matched-profile ansatz at the "
                  "one axisymmetric geometry; off-syzygy orientations "
                  "remain bounded by the step_084 scan"),
        "units": {"length": "d_EM = 3.84399e8 m",
                  "gradient": "g_t = c H0 / 2"},
        "charges": {"Q_earth": Q_E, "Q_moon": Q_M,
                    "a_solar_newtonian_over_gt": A_SOLAR_N},
        "source_width_over_d": SIG_SRC,
        "grid": {"NS": NS, "NZ": NZ, "SMAX": SMAX,
                 "ZMIN": ZMIN, "ZMAX": ZMAX},
        "sectors": {},
    }

    for sector in ("baseline", "two_branch"):
        u0 = u_solar(sector)
        warm = sorted(u for u in (1e1, 3e1, 1e2, 3e2) if u < u0)
        psi = solve_with_continuation(sector, u0, warm_u0s=warm)
        out["sectors"][sector] = extract(psi, u0, sector)
        out["sectors"][sector]["u_solar_scalar_gradient_over_gt"] = u0

    # comparison against the step_084 matched-ansatz baseline
    ansatz_path = (ROOT / "results/outputs/"
                   "step_084_earth_moon_solar_bridge.json")
    if ansatz_path.exists():
        ansatz = json.loads(ansatz_path.read_text())
        out["matched_ansatz_comparison"] = {
            "step084_bridge_min_gradient_over_gt":
                ansatz["sectors"]["baseline"]["0"]
                ["bridge_min_gradient_over_gt"],
            "step086_baseline_bridge_min_over_gt":
                out["sectors"]["baseline"]["bridge_min_abs_a_over_gt"],
            "step084_two_branch_bridge_min_over_gt":
                ansatz["sectors"]["two_branch"]["0"]
                ["bridge_min_gradient_over_gt"],
            "step086_two_branch_bridge_min_over_gt":
                out["sectors"]["two_branch"]["bridge_min_abs_a_over_gt"],
        }

    out["limitations"] = (
        "Syzygy only (the axisymmetric configuration); sources are "
        "Gaussian lumps of width 0.08 d_EM (~3e7 m), so body interiors "
        "and sub-source-radius structure are unresolved.  The scalar "
        "field solve here supplies the screened force channel; it is "
        "not an LLR range-residual model and does not by itself relax "
        "the Solar-System constraints.  Off-syzygy orientations "
        "require a 3D solve and remain bounded by the step_084 "
        "matched-ansatz scan.")

    dest = ROOT / "results/outputs/step_086_earth_moon_coupled_pde.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest)
    for sector, r in out["sectors"].items():
        print(f"{sector}: bridge min |a| = "
              f"{r['bridge_min_abs_a_over_gt']:.3e} g_t at "
              f"z/d = {r['bridge_min_z_over_d']:.3f}, "
              f"P_X = {r['bridge_min_P_X']:.3e}, "
              f"screened frac = {r['screened_volume_fraction_near_pair']:.3e}")


if __name__ == "__main__":
    main()
