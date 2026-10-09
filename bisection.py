# Reference: Numerical Methods for Eigenvalue Problems (2012) Ch. 6
# Additional reference for sturm sequences: https://www.cse.psu.edu/~b58/cse456/lecture13.pdf

# matrix must be real symmetric tridiagonal

import numpy as np


def sturm_evaluate(z, d, e):
    """Count sign changes in q_i(z) = det(z I - T_i) for entries d, e."""
    q = [1, z - d[0]]
    count = 0
    if q[0] * q[1] < 0 or q[1] == 0:
        count += 1
    n = len(d)
    for i in range(2, n + 1):
        q.append((z - d[i - 1]) * q[i - 1] - abs(e[i - 2]) ** 2 * q[i - 2])
        if q[i] * q[i - 1] < 0 or q[i] == 0:
            count += 1
    return count


def sturm_bisection(index, d, e, lower, upper, tol=1e-6):
    """Approximate the one-based index-th eigenvalue of tridiagonal T."""
    n = d.shape[0]
    while upper - lower > tol:
        midpoint = (upper + lower) / 2
        count = sturm_evaluate(midpoint, d, e)
        print(index, count)

        if index <= n - count:
            upper = midpoint
        else:
            lower = midpoint
    print()
    return midpoint


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
    T = symmetric_tridiagonal_matrix(10)
    d = np.diag(T)
    e = np.diag(T, -1)

    lower, upper = gershgorin_bound(d, e)

    n = T.shape[0]
    eigvals = []

    for index in range(1, n + 1):
        eigval = sturm_bisection(index, d, e, lower, upper)
        eigvals.append(eigval)
        print()

    print(eigvals)
    print(sorted(np.linalg.eig(T)[0]))


if __name__ == "__main__":
    main()
