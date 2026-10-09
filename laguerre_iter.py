"""Laguerre root finding for characteristic polynomials of small matrices."""

import numpy as np

from _validation import iteration_limit, tolerance, tridiagonal_entries
from bisection import gershgorin_bound, symmetric_tridiagonal_matrix


def construct_characteristic_polynomial(d, e):
    """Return p_n(z) = det(T - z I) for real symmetric tridiagonal T.

    p_0 = 1, p_1 = d_1 - z, and
    p_i(z) = (d_i - z) p_{i-1}(z) - e_{i-1}^2 p_{i-2}(z).
    """
    d, e = tridiagonal_entries(d, e)
    p_prev, p = np.poly1d([1.]), np.poly1d([-1., d[0]])
    for i in range(1, len(d)):
        p_next = np.polymul(p, [-1., d[i]]) - e[i - 1]**2 * p_prev
        p_prev, p = p, p_next
    if not np.all(np.isfinite(p.c)):
        raise FloatingPointError("Characteristic polynomial coefficients overflowed")
    return p


def laguerre_method(p, z0, tol=1e-6, max_iter=100):
    """Find one real or complex root of a poly1d or coefficient array.

    Use the polynomial degree and the larger-magnitude Laguerre denominator.
    Convergence uses a coefficient-scaled polynomial residual; exhausted
    iterations raise RuntimeError. High-degree polynomial deflation can be
    ill-conditioned even when an individual polynomial root has converged.
    """
    tol, max_iter = tolerance(tol), iteration_limit(max_iter)
    p = np.poly1d(p)
    if len(p) < 1 or not np.all(np.isfinite(p.c)):
        raise ValueError("p must be a finite nonconstant polynomial")
    if not np.isscalar(z0) or not np.isfinite(z0):
        raise ValueError("z0 must be a finite scalar")
    p = np.poly1d(p.c / np.max(np.abs(p.c)))
    degree = len(p)
    p_prime, p_double_prime = np.polyder(p), np.polyder(p, 2)
    z = complex(z0)
    for k in range(max_iter + 1):
        p_value = np.polyval(p, z)
        evaluation_bound = np.polyval(np.abs(p.c), abs(z))
        if not np.isfinite(p_value) or not np.isfinite(evaluation_bound):
            raise FloatingPointError("Polynomial evaluation exceeded floating-point range")
        if abs(p_value) <= tol * evaluation_bound:
            return np.real_if_close(z).item()
        if k == max_iter:
            break
        log_derivative = np.polyval(p_prime, z) / p_value
        log_curvature = log_derivative**2 - np.polyval(p_double_prime, z) / p_value
        radical = np.emath.sqrt((degree - 1) * (degree * log_curvature - log_derivative**2))
        plus, minus = log_derivative + radical, log_derivative - radical
        denominator = plus if abs(plus) >= abs(minus) else minus
        if denominator == 0:
            # Escape a stationary point without an undefined division.
            step = (1 + abs(z)) * np.exp(1j * (k + 1))
        else:
            step = degree / denominator
        if not np.isfinite(step):
            raise FloatingPointError("Laguerre iteration produced an invalid step")
        z -= step
    raise RuntimeError("Laguerre iteration did not converge within max_iter")


def main():
    # Coefficient formation and deflation are intended for small examples.
    d = np.array([-3., -1., 1., 2., 4.])
    e = np.array([.2, -.3, .1, .4])
    T = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
    lower, upper = gershgorin_bound(d, e)
    p = construct_characteristic_polynomial(d, e)
    eigvals = []
    for index in range(len(d)):
        eigval = laguerre_method(p, (lower + upper) / 2, tol=1e-12)
        eigvals.append(eigval)
        p = np.polydiv(p, [-1., eigval])[0]
    print(np.sort(np.real_if_close(eigvals)))
    print(np.linalg.eigvalsh(T))


if __name__ == "__main__":
    main()
