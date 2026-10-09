# Reference: Numerical Methods for Eigenvalue Problems (2012) Ch. 6

# matrix must be real symmetric tridiagonal

import numpy as np


def construct_characteristic_polynomial(d, e):
    """Return p_n(z) = det(T - z I) for real symmetric tridiagonal T.

    d contains the n diagonal entries; e contains the n - 1 off-diagonal
    entries. For leading principal blocks T_i, the recurrence is
    p_i(z) = (d_i - z) p_{i-1}(z) - e_{i-1}^2 p_{i-2}(z).
    """
    n = d.shape[0]
    p = [1, np.poly1d([-1, d[0]])]
    for i in range(1, n):
        diagonal_term = np.polymul(p[i], [-1, d[i]])
        offdiag_term = np.polymul(p[i - 1], [-e[i - 1]**2])
        p.append(np.polyadd(diagonal_term, offdiag_term))
    return p[-1]


def laguerre_method(p, z0, tol=1e-6, max_iter=100):
    """
    Implements the Laguerre method for finding a root of a polynomial function.

    Parameters:
    p (np.poly1d): Polynomial whose root is sought.
    z0 (float): Initial guess for the root.
    tol (float): Tolerance for convergence.
    max_iter (int): Maximum number of iterations.

    Returns:
    float: Approximation of the root.
    """
    z = z0
    degree = len(p)

    p_prime = np.polyder(p)
    p_double_prime = np.polyder(p_prime)

    for k in range(max_iter):
        p_value = np.polyval(p, z)
        p_prime_value = np.polyval(p_prime, z)
        p_double_prime_value = np.polyval(p_double_prime, z)

        if abs(p_value) < tol:
            return z

        log_derivative = p_prime_value / p_value
        log_curvature = log_derivative**2 - p_double_prime_value / p_value

        if log_derivative >= 0:
            step = degree / (
                log_derivative
                + np.emath.sqrt((degree - 1) * (degree * log_curvature - log_derivative**2))
            )
        else:
            step = degree / (
                log_derivative
                - np.emath.sqrt((degree - 1) * (degree * log_curvature - log_derivative**2))
            )

        z -= step

        if abs(step) < tol:
            return z

    Exception('Maximum number of iterations exceeded.')


def gershgorin_bound(d, e):
    """Return lower and upper eigenvalue bounds for tridiagonal entries d, e."""
    n = d.shape[0]
    lower = np.min(
        [d[0] - np.abs(e[0]), d[n - 1] - np.abs(e[n - 2])]
        + [d[i] - np.abs(e[i]) - np.abs(e[i - 1]) for i in range(1, n - 1)]
    )
    upper = np.max(
        [d[0] + np.abs(e[0]), d[n - 1] + np.abs(e[n - 2])]
        + [d[i] + np.abs(e[i]) + np.abs(e[i - 1]) for i in range(1, n - 1)]
    )
    return lower, upper


def symmetric_tridiagonal_matrix(n):
    """Generate an n x n real symmetric tridiagonal matrix T."""
    d = np.random.rand(n)
    e = np.random.rand(n-1)
    T = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
    return T


def main():
    T = symmetric_tridiagonal_matrix(30)
    d = np.diag(T)
    e = np.diag(T, -1)

    lower, upper = gershgorin_bound(d, e)

    n = T.shape[0]
    eigvals = []

    p = construct_characteristic_polynomial(d, e)
    z0 = (lower + upper) / 2

    for index in range(1, n + 1):
        eigval = laguerre_method(p, z0)
        eigvals.append(eigval)

        # update the characteristic polynomial
        p = np.polydiv(p, [-1, eigval])[0]

    print(sorted(eigvals))
    print(sorted(np.linalg.eig(T)[0]))


if __name__ == "__main__":
    main()
