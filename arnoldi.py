import numpy as np


def arnoldi_iteration(A, x0, m: int, tol=1e-12):
    """Compute a basis of span{x0, A x0, ..., A^m x0} for real A.

    Arguments
      A: n x n array
      x0: initial vector (length n)
      m: number of Arnoldi steps, must be >= 1
      tol: tolerance for detecting Krylov breakdown

    Returns
      Q: n x (m + 1) array containing the basis vectors
      H: (m + 1) x m upper Hessenberg representation of A

    On early breakdown, unused columns of Q and entries of H remain zero.
    """
    n = A.shape[0]
    H = np.zeros((m + 1, m))
    Q = np.zeros((n, m + 1))
    # Normalize the input vector
    Q[:, 0] = x0 / np.linalg.norm(x0, 2)
    for k in range(1, m + 1):
        y = np.dot(A, Q[:, k - 1])  # Generate a new candidate vector
        for j in range(k):  # Subtract the projections on previous vectors
            H[j, k - 1] = np.dot(Q[:, j].T, y)
            y = y - H[j, k - 1] * Q[:, j]
        H[k, k - 1] = np.linalg.norm(y, 2)
        if H[k, k - 1] > tol:
            Q[:, k] = y / H[k, k - 1]
        else:  # Stop when the next basis vector is too small to normalize.
            return Q, H
    return Q, H
