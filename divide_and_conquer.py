import numpy as np

from _validation import tridiagonal_matrix
from bisection import symmetric_tridiagonal_matrix

# Reference: https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf


def partition(T):
    """Return block diagonal B and its rank-one update without modifying T."""
    T = tridiagonal_matrix(T)
    n = len(T)
    if n < 2:
        raise ValueError("Partition requires at least two rows")
    split = n // 2
    rho = T[split - 1, split]
    B = T.copy()
    B[split - 1, split - 1] -= rho
    B[split, split] -= rho
    B[split - 1, split] = B[split, split - 1] = 0
    v = np.zeros(n)
    v[split - 1:split + 1] = 1
    return B, rho * np.outer(v, v)


def div_conq(T):
    """Return sorted eigenvalues and eigenvector columns Q of real tridiagonal T.

    The rank-one merge uses NumPy's dense symmetric eigensolver. All calls
    preserve their input and return a one-dimensional eigenvalue array.
    """
    T = tridiagonal_matrix(T)
    n = len(T)
    if n == 1:
        return np.array([T[0, 0]]), np.eye(1)
    split = n // 2
    B, _ = partition(T)
    eigvals1, Q1 = div_conq(B[:split, :split])
    eigvals2, Q2 = div_conq(B[split:, split:])
    D = np.diag(np.concatenate((eigvals1, eigvals2)))
    v = np.concatenate((Q1[-1], Q2[0]))
    rho = T[split - 1, split]
    eigvals, Q_update = np.linalg.eigh(D + rho * np.outer(v, v))
    Q = np.zeros((n, n))
    Q[:split, :split] = Q1
    Q[split:, split:] = Q2
    return eigvals, Q @ Q_update


def main():
    T = np.diag([2., 1., 3., 4., 5.])
    e = np.array([.5, -.5, .25, .5])
    T += np.diag(e, 1) + np.diag(e, -1)
    eigvals, Q = div_conq(T)
    print(eigvals)
    print(np.linalg.eigvalsh(T))


if __name__ == "__main__":
    main()
