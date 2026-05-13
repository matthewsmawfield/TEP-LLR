"""
Stable linear-algebra defaults for LLR-scale design matrices.

Julian-date harmonics and stacked periodic controls can be nearly rank-deficient;
``numpy.linalg.lstsq(..., rcond=None)`` can trigger noisy LAPACK matmul overflows
on some platforms.  Pipeline code should call ``stable_lstsq`` instead.
"""

from __future__ import annotations

import contextlib
import warnings

import numpy as np

# Relative cutoff for singular values (max(s) * rcond); stable for n~26k, k~30.
_LSTSQ_RCOND = 1e-8


@contextlib.contextmanager
def suppress_scipy_array_api_matmul_runtime_warning():
    """
    SciPy 1.14+ can emit benign RuntimeWarning from array-api matmul paths
    inside ``pearsonr`` / ``corrcoef`` on large float64 vectors (NumPy 2.x).
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=RuntimeWarning,
            message=".*encountered in matmul.*",
        )
        yield


def stable_lstsq(
    a: np.ndarray,
    b: np.ndarray,
    *,
    rcond: float | None = None,
) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:
    """
    Ordinary least squares with float64 inputs and a finite singular-value cutoff.

    Parameters match ``numpy.linalg.lstsq`` return tuple for drop-in replacement.
    """
    a = np.asarray(a, dtype=np.float64, order="C")
    b = np.asarray(b, dtype=np.float64, order="C")
    if b.ndim == 1:
        b = b.reshape(-1)
    rc = _LSTSQ_RCOND if rcond is None else float(rcond)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        return np.linalg.lstsq(a, b, rcond=rc)


def hat_diagonal_from_qr(X: np.ndarray) -> np.ndarray:
    """
    Leverage h_ii = (H)_ii for H = X (X'X)^{-1} X' without forming H.

    Uses reduced QR: for X = Q R with Q n×k orthonormal columns,
    H = Q Q' and h_ii = sum_j Q[i,j]^2.
    """
    X = np.asarray(X, dtype=np.float64, order="C")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        Q, _ = np.linalg.qr(X, mode="reduced")
    return np.sum(Q * Q, axis=1)
