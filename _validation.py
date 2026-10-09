"""Shared input checks for the numerical examples."""

from numbers import Integral

import numpy as np


def matrix(A, *, square=True, symmetric=False, real=False):
    A = np.asarray(A)
    if A.ndim != 2 or 0 in A.shape:
        raise ValueError("A must be a nonempty two-dimensional matrix")
    if square and A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if real and np.iscomplexobj(A):
        if np.any(A.imag != 0):
            raise ValueError("This routine requires a real matrix")
        A = A.real
    A = np.array(A, dtype=complex if np.iscomplexobj(A) else float, copy=True)
    if not np.all(np.isfinite(A)):
        raise ValueError("A must contain finite entries")
    if symmetric:
        scale = max(np.max(np.abs(A)), np.finfo(float).tiny)
        if not np.allclose(A / scale, A.conj().T / scale, rtol=1e-12, atol=1e-14):
            raise ValueError("A must be symmetric or Hermitian")
    return A


def vector(x, n, *, nonzero=True):
    x = np.asarray(x)
    if x.shape != (n,) or not np.all(np.isfinite(x)):
        raise ValueError(f"The vector must have {n} finite entries")
    x = np.array(x, dtype=complex if np.iscomplexobj(x) else float, copy=True)
    if nonzero and np.linalg.norm(x) == 0:
        raise ValueError("The initial vector must be nonzero")
    return x


def tolerance(tol):
    if not np.isscalar(tol) or np.iscomplexobj(tol) or not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be a positive finite number")
    return float(tol)


def iteration_limit(max_iter):
    if isinstance(max_iter, bool) or not isinstance(max_iter, Integral) or max_iter < 0:
        raise ValueError("max_iter must be a nonnegative integer")
    return int(max_iter)


def spectral_shift(shift):
    if not np.isscalar(shift) or not np.isfinite(shift):
        raise ValueError("shift must be a finite scalar")
    return shift


def tridiagonal_entries(d, e):
    d, e = np.asarray(d), np.asarray(e)
    if d.ndim != 1 or len(d) == 0 or e.shape != (len(d) - 1,):
        raise ValueError("d must have n entries and e must have n - 1 entries")
    if np.iscomplexobj(d) or np.iscomplexobj(e):
        raise ValueError("Tridiagonal entries must be real")
    d, e = np.array(d, dtype=float, copy=True), np.array(e, dtype=float, copy=True)
    if not np.all(np.isfinite(d)) or not np.all(np.isfinite(e)):
        raise ValueError("Tridiagonal entries must be finite")
    return d, e


def tridiagonal_matrix(T):
    T = matrix(T, symmetric=True, real=True)
    if np.any(np.triu(T, 2) != 0) or np.any(np.tril(T, -2) != 0):
        raise ValueError("T must be tridiagonal")
    return T
