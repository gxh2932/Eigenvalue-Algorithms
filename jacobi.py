import numpy as np

from _validation import iteration_limit, matrix, tolerance


def jacobi_eigenvalue_algorithm(A, tol=1e-10, max_iter=10000):
    """Return eigvals, eigenvector columns Q, and n_iter for real symmetric A.

    tol is relative to the matrix norm. Raise RuntimeError on exhaustion.
    """
    A = matrix(A, symmetric=True, real=True)
    tol, max_iter = tolerance(tol), iteration_limit(max_iter)
    n = len(A)
    Q = np.eye(n)
    scale = np.linalg.norm(A, ord=np.inf)
    if scale == 0:
        return np.zeros(n), Q, 0
    A /= scale
    for n_iter in range(max_iter + 1):
        offdiag = np.abs(np.triu(A, 1))
        i, j = np.unravel_index(np.argmax(offdiag), A.shape)
        if offdiag[i, j] <= tol:
            return np.diag(A) * scale, Q, n_iter
        if n_iter == max_iter:
            break
        theta = .5 * np.arctan2(2 * A[i, j], A[i, i] - A[j, j])
        c, s = np.cos(theta), np.sin(theta)
        J = np.eye(n)
        J[i, i] = J[j, j] = c
        J[i, j], J[j, i] = -s, s
        A = J.T @ A @ J
        A = (A + A.T) / 2
        Q = Q @ J
    raise RuntimeError("Jacobi iteration did not converge within max_iter")


def main():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(5, 5))
    A = A + A.T
    eigvals, Q, n_iter = jacobi_eigenvalue_algorithm(A)
    print(np.sort(eigvals))
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
