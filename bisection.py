"""Sturm bisection for real symmetric tridiagonal matrices."""

import numpy as np

from _validation import iteration_limit, tolerance, tridiagonal_entries


def sturm_evaluate(z, d, e):
    """Return the number of eigenvalues strictly greater than z.

    Use LDL pivots of T - z I, rather than products of characteristic
    polynomials. A zero pivot is replaced by a small negative pivot, so an
    eigenvalue equal to z belongs to the lower (<= z) part of the spectrum.
    Zero off-diagonal entries naturally start a new independent block.
    """
    d, e = tridiagonal_entries(d, e)
    if not np.isscalar(z) or np.iscomplexobj(z) or np.isnan(z):
        raise ValueError("z must be a real scalar")
    if np.isposinf(z):
        return 0
    if np.isneginf(z):
        return len(d)
    scale = max(np.max(np.abs(d)), np.max(np.abs(e), initial=0), abs(z))
    if scale == 0:
        return 0
    # Subtract close, same-sign numbers before scaling to retain their gap.
    # For opposite signs, scale first to avoid overflowing d - z.
    same_sign = np.signbit(d) == np.signbit(z)
    shifted_diagonal = np.empty_like(d)
    shifted_diagonal[same_sign] = (d[same_sign] - z) / scale
    shifted_diagonal[~same_sign] = d[~same_sign] / scale - z / scale
    e = e / scale
    pivmin = np.finfo(float).tiny
    pivot = shifted_diagonal[0]
    count_le = 0
    for i in range(len(d)):
        if i:
            pivot = shifted_diagonal[i] - (e[i - 1] ** 2) / pivot
        if abs(pivot) < pivmin:
            pivot = -pivmin
        count_le += int(pivot < 0)
    return len(d) - count_le


def sturm_bisection(index, d, e, lower, upper, tol=1e-6, max_iter=1000):
    """Approximate the one-based index-th eigenvalue; tol is absolute."""
    d, e = tridiagonal_entries(d, e)
    tol, max_iter = tolerance(tol), iteration_limit(max_iter)
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)) or not 1 <= index <= len(d):
        raise ValueError("index must be an integer between 1 and n")
    if not np.isfinite(lower) or not np.isfinite(upper) or lower > upper:
        raise ValueError("lower and upper must be finite ordered bounds")
    if (len(d) - sturm_evaluate(np.nextafter(lower, -np.inf), d, e) >= index
            or len(d) - sturm_evaluate(upper, d, e) < index):
        raise ValueError("The interval does not bracket the requested eigenvalue")
    for k in range(max_iter + 1):
        midpoint = lower / 2 + upper / 2
        if upper - lower <= tol or midpoint == lower or midpoint == upper:
            return midpoint
        if k == max_iter:
            break
        count_le = len(d) - sturm_evaluate(midpoint, d, e)
        if index <= count_le:
            upper = midpoint
        else:
            lower = midpoint
    raise RuntimeError("Bisection did not converge within max_iter")


def gershgorin_bound(d, e):
    """Return enclosing eigenvalue bounds, including for a 1 x 1 matrix."""
    d, e = tridiagonal_entries(d, e)
    radii = np.zeros(len(d))
    radii[:-1] += np.abs(e)
    radii[1:] += np.abs(e)
    lower, upper = np.min(d - radii), np.max(d + radii)
    if not np.isfinite(lower) or not np.isfinite(upper):
        raise FloatingPointError("Eigenvalue bounds exceed floating-point range")
    return lower, upper


def symmetric_tridiagonal_matrix(n):
    """Generate an n x n real symmetric tridiagonal matrix T."""
    n = iteration_limit(n)
    if n == 0:
        raise ValueError("n must be positive")
    d = np.random.rand(n)
    e = np.random.rand(n - 1)
    return np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)


def main():
    d = np.array([2., 1., 3.])
    e = np.zeros(2)
    T = np.diag(d)
    lower, upper = gershgorin_bound(d, e)
    eigvals = [sturm_bisection(i, d, e, lower, upper) for i in range(1, len(d) + 1)]
    print(eigvals)
    print(np.linalg.eigvalsh(T))


if __name__ == "__main__":
    main()
