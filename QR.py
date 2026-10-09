import numpy as np


def QR_decomposition(A):
    """Return the orthonormal basis Q and upper triangular factor R of A."""
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    for j in range(n):
        y = A[:, j]
        for i in range(j):
            R[i, j] = np.dot(Q[:, i], A[:, j])
            y = y - R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(y)
        Q[:, j] = y / R[j, j]

    return Q, R


def QR(A, max_iter=1000):
    """Approximate the eigenvalues of real symmetric A using QR iteration."""
    for k in range(max_iter):
        Q, R = QR_decomposition(A)
        A = R @ Q

    return np.diag(A)


def main():
    A = np.random.rand(3, 3)
    A = A + A.T

    print(np.sort(QR(A)))
    print(np.sort(np.linalg.eigvals(A)))


if __name__ == "__main__":
    main()
