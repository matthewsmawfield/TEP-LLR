#!/usr/bin/env python3
"""Step 077: Formal GP vs cosD model comparison on held-out data.

Adjudicates between the synodic cosD model and the 2D GP on elongation x time
using out-of-sample predictive performance (cross-validated log-likelihood)
rather than in-sample AIC alone.
"""

from pathlib import Path
import json
import sys
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.model_selection import KFold

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "utils"))

from logger import TEPLogger

def load_step_061():
    path = PROJECT_ROOT / "results" / "outputs" / "step_061_systematic_sensitivity_analysis.json"
    with open(path, "r") as f:
        return json.load(f)

def build_gp_data():
    """Reconstruct the 2D binned elongation x time data from step_061 output."""
    step_061 = load_step_061()
    gp_info = step_061["adversarial_gp"]
    
    # Load processed data to reconstruct bins
    import pandas as pd
    df = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv")
    
    # Bin by elongation (12 bins) and time (yearly bins)
    # elongation is in radians; convert to degrees for 12 bins
    df["elong_deg"] = np.degrees(df["elongation_rad"])
    df["elong_bin"] = pd.cut(df["elong_deg"], bins=12, labels=False)
    df["year_bin"] = pd.cut(df["date_julian_year"], bins=12, labels=False)
    
    # Compute bin means
    bin_means = df.groupby(["elong_bin", "year_bin"])["residual_m"].mean().reset_index()
    bin_means = bin_means.dropna()
    
    X = bin_means[["elong_bin", "year_bin"]].values.astype(float)
    y = bin_means["residual_m"].values
    
    return X, y, gp_info

def cosd_model_loglik(X, y, elong_bins=12):
    """Fit cosD model on binned data and return log-likelihood."""
    # elong_bin corresponds to phase; cosD = cos(2*pi * elong_bin / elong_bins)
    cosD = np.cos(2 * np.pi * X[:, 0] / elong_bins)
    # OLS fit
    A = np.vstack([cosD, np.ones(len(cosD))]).T
    beta, residuals, rank, s = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ beta
    rss = np.sum((y - y_pred) ** 2)
    n = len(y)
    sigma2 = rss / n
    loglik = -0.5 * n * np.log(2 * np.pi * sigma2) - 0.5 * rss / sigma2
    return loglik, beta, sigma2

def gp_model_loglik(X_train, y_train, X_test, y_test, kernel):
    """Fit GP on train, return log-likelihood on test."""
    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=3)
    gp.fit(X_train, y_train)
    y_pred, y_std = gp.predict(X_test, return_std=True)
    # Gaussian log-likelihood on test data
    loglik = -0.5 * np.sum(((y_test - y_pred) / y_std) ** 2) - np.sum(np.log(y_std)) - 0.5 * len(y_test) * np.log(2 * np.pi)
    return loglik, gp

def run_gp_model_comparison():
    logger = TEPLogger("step_077_gp_formal_model_comparison")
    
    X, y, gp_info = build_gp_data()
    n = len(X)
    
    # Kernel from step_061
    kernel = ConstantKernel(0.1**2) * RBF(length_scale=[1.74, 1.09]) + WhiteKernel(noise_level=0.000107)
    
    # 5-fold CV for robustness
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    
    cosd_lls = []
    gp_lls = []
    joint_lls = []
    
    for train_idx, test_idx in kf.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # cosD-only model
        cosD_train = np.cos(2 * np.pi * X_train[:, 0] / 12)
        A_train = np.vstack([cosD_train, np.ones(len(cosD_train))]).T
        beta, _, _, _ = np.linalg.lstsq(A_train, y_train, rcond=None)
        
        cosD_test = np.cos(2 * np.pi * X_test[:, 0] / 12)
        A_test = np.vstack([cosD_test, np.ones(len(cosD_test))]).T
        y_pred_cosd = A_test @ beta
        sigma2 = np.mean((y_train - A_train @ beta) ** 2)
        ll_cosd = -0.5 * np.sum(((y_test - y_pred_cosd) / np.sqrt(sigma2)) ** 2) - 0.5 * len(y_test) * np.log(2 * np.pi * sigma2)
        cosd_lls.append(ll_cosd)
        
        # GP model
        ll_gp, gp = gp_model_loglik(X_train, y_train, X_test, y_test, kernel)
        gp_lls.append(ll_gp)
        
        # Joint model: cosD + GP mean as additional predictor
        gp.fit(X_train, y_train)
        gp_mean_train = gp.predict(X_train)
        A_joint_train = np.vstack([cosD_train, gp_mean_train, np.ones(len(cosD_train))]).T
        beta_joint, _, _, _ = np.linalg.lstsq(A_joint_train, y_train, rcond=None)
        
        gp_mean_test = gp.predict(X_test)
        A_joint_test = np.vstack([cosD_test, gp_mean_test, np.ones(len(cosD_test))]).T
        y_pred_joint = A_joint_test @ beta_joint
        sigma2_joint = np.mean((y_train - A_joint_train @ beta_joint) ** 2)
        ll_joint = -0.5 * np.sum(((y_test - y_pred_joint) / np.sqrt(sigma2_joint)) ** 2) - 0.5 * len(y_test) * np.log(2 * np.pi * sigma2_joint)
        joint_lls.append(ll_joint)
    
    result = {
        "step_id": "step_077",
        "status": "PASS",
        "n_bins": n,
        "n_folds": 5,
        "models": {
            "cosd_only": {
                "mean_cv_loglik": float(np.mean(cosd_lls)),
                "std_cv_loglik": float(np.std(cosd_lls)),
                "logliks": [float(x) for x in cosd_lls]
            },
            "gp_only": {
                "mean_cv_loglik": float(np.mean(gp_lls)),
                "std_cv_loglik": float(np.std(gp_lls)),
                "logliks": [float(x) for x in gp_lls]
            },
            "joint_cosd_plus_gp": {
                "mean_cv_loglik": float(np.mean(joint_lls)),
                "std_cv_loglik": float(np.std(joint_lls)),
                "logliks": [float(x) for x in joint_lls]
            }
        },
        "comparison": {
            "gp_vs_cosd_delta_ll": float(np.mean(gp_lls) - np.mean(cosd_lls)),
            "joint_vs_cosd_delta_ll": float(np.mean(joint_lls) - np.mean(cosd_lls)),
            "joint_vs_gp_delta_ll": float(np.mean(joint_lls) - np.mean(gp_lls)),
            "gp_wins_cv": bool(np.mean(gp_lls) > np.mean(cosd_lls)),
            "preferred_model_cv": "gp_only" if np.mean(gp_lls) > np.mean(joint_lls) and np.mean(gp_lls) > np.mean(cosd_lls) else (
                "joint_cosd_plus_gp" if np.mean(joint_lls) > np.mean(cosd_lls) and np.mean(joint_lls) > np.mean(gp_lls) else "cosd_only"
            )
        },
        "interpretation": "5-fold CV predictive log-likelihood comparison. If GP-only wins, the non-parametric surface generalises better than the synodic sinusoid on held-out data. If cosD-only or joint wins, the synodic component has predictive validity."
    }
    
    output_path = PROJECT_ROOT / "results" / "outputs" / "step_077_gp_formal_model_comparison.json"
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)
    
    logger.info(f"Step 077 complete. Preferred model (CV): {result['comparison']['preferred_model_cv']}")
    logger.info(f"GP vs cosD delta LL: {result['comparison']['gp_vs_cosd_delta_ll']:.2f}")
    
    return result

if __name__ == "__main__":
    run_gp_model_comparison()
