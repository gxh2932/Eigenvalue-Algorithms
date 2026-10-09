import numpy as np

from _validation import iteration_limit, matrix, tolerance


def QR_decomposition(A):
    """Return reduced Householder QR factors, including for rank-deficient A."""
    A = matrix(A, square=False)
    return np.linalg.qr(A, mode="reduced")


def QR(A, max_iter=1000, tol=1e-12):
    """Return eigenvalues of real symmetric A by shifted QR with deflation.

    tol is relative to the matrix norm. Raise RuntimeError on nonconvergence.
    """
    A = matrix(A, symmetric=True, real=True)
    max_iter, tol = iteration_limit(max_iter), tolerance(tol)
    n = A.shape[0]
    eigvals = np.empty(n)
    scale = np.linalg.norm(A, ord=np.inf)
    if scale == 0:
        return np.zeros(n)
    A /= scale
    active, n_iter = n, 0
    while active > 1:
        if np.linalg.norm(A[active - 1, :active - 1]) <= tol:
            eigvals[active - 1] = A[active - 1, active - 1]
            active -= 1
            continue
        if n_iter >= max_iter:
            raise RuntimeError("QR iteration did not converge within max_iter")
        # Wilkinson shift from the trailing 2 x 2 principal block.
        a, b = A[active - 2, active - 2], A[active - 2, active - 1]
        c = A[active - 1, active - 1]
        delta = (a - c) / 2
        sign = 1.0 if delta >= 0 else -1.0
        denominator = abs(delta) + np.hypot(delta, b)
        shift = c if denominator == 0 else c - sign * b * (b / denominator)
        Q, R = QR_decomposition(A[:active, :active] - shift * np.eye(active))
        block = R @ Q + shift * np.eye(active)
        A[:active, :active] = (block + block.T) / 2
        n_iter += 1
    eigvals[0] = A[0, 0]
    return eigvals * scale


def main():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(5, 5))
    A = A + A.T
    print(np.sort(QR(A)))
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
