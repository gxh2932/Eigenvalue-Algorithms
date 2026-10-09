import numpy as np

from _validation import iteration_limit, matrix, tolerance, vector


def tridiag(e_lower, d, e_upper, lower_offset=-1, diagonal_offset=0, upper_offset=1):
    """Construct T from its lower, main, and upper diagonal entries."""
    return (np.diag(e_lower, lower_offset) + np.diag(d, diagonal_offset)
            + np.diag(e_upper, upper_offset))


def lanczos(A, m=None, tol=1e-12, x0=None):
    """Return a real tridiagonal Ritz matrix for symmetric or Hermitian A.

    Run at most m steps (default n), reorthogonalizing the basis. Stop before
    normalizing a zero residual. On breakdown, T is smaller than n x n and
    describes the invariant subspace reached from x0, not all multiplicities.
    tol is relative to ||A||.
    """
    A = matrix(A, symmetric=True)
    n = len(A)
    m = n if m is None else iteration_limit(m)
    tol = tolerance(tol)
    if not 1 <= m <= n:
        raise ValueError("m must be between 1 and n")
    x0 = np.ones(n) if x0 is None else vector(x0, n)
    q = x0 / np.linalg.norm(x0)
    q_prev = np.zeros(n)
    beta = 0.
    d, e, basis = [], [], []
    scale = max(np.linalg.norm(A, ord=np.inf), np.finfo(float).tiny)
    for k in range(m):
        basis.append(q.copy())
        y = A @ q
        alpha = np.vdot(q, y).real
        r = y - alpha * q - beta * q_prev
        Q = np.column_stack(basis)
        for _ in range(2):
            r -= Q @ (Q.conj().T @ r)
        beta = np.linalg.norm(r)
        d.append(alpha)  # Preserve the signed Rayleigh coefficient.
        if beta <= tol * scale or k == m - 1:
            break
        e.append(beta)
        q_prev, q = q, r / beta
    return tridiag(e, d, e)


def generate_hermite_matrix(n):
    """Generate the n x n real symmetric matrix A used in the example."""
    n = iteration_limit(n)
    if n == 0:
        raise ValueError("n must be positive")
    return np.diag(2. * np.arange(n)) + np.diag(-np.ones(n - 1), 1) + np.diag(-np.ones(n - 1), -1)


def main():
    A = np.diag([-3., -1., 2.])
    T = lanczos(A)
    print(np.linalg.eigvalsh(T))
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
