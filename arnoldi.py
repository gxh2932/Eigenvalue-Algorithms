import numpy as np

from _validation import iteration_limit, matrix, tolerance, vector


def arnoldi_iteration(A, x0, m: int, tol=1e-12):
    """Return Q and H satisfying A Q[:, :m] = Q H for a Krylov basis.

    A is n x n; m is the number of steps (1 <= m <= n). Q has at most
    m + 1 columns, and H has one column per completed step. On breakdown,
    return only the active basis and square H. tol is relative to ||A||.
    """
    A = matrix(A)
    m, tol = iteration_limit(m), tolerance(tol)
    n = len(A)
    if not 1 <= m <= n:
        raise ValueError("m must be between 1 and n")
    x0 = vector(x0, n)
    dtype = np.result_type(A, x0)
    Q = np.zeros((n, m + 1), dtype=dtype)
    H = np.zeros((m + 1, m), dtype=dtype)
    Q[:, 0] = x0 / np.linalg.norm(x0)
    scale = max(np.linalg.norm(A, ord=np.inf), np.finfo(float).tiny)
    for k in range(m):
        y = A @ Q[:, k]
        # Twice-modified Gram-Schmidt keeps the computed basis orthogonal.
        for _ in range(2):
            coefficients = Q[:, :k + 1].conj().T @ y
            H[:k + 1, k] += coefficients
            y -= Q[:, :k + 1] @ coefficients
        beta = np.linalg.norm(y)
        if beta <= tol * scale or k + 1 == n:
            return Q[:, :k + 1], H[:k + 1, :k + 1]
        H[k + 1, k] = beta
        Q[:, k + 1] = y / beta
    return Q, H
