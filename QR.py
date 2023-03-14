import numpy as np


# real symmetric matrix


def QR_decomposition(A):
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    for k in range(n):
        v = A[:, k]
        for j in range(k):
            R[j, k] = np.dot(Q[:, j], A[:, k])
            v = v - R[j, k] * Q[:, j]
        R[k, k] = np.linalg.norm(v)
        Q[:, k] = v / R[k, k]

    return Q, R


def QR(A, num_iters=1000):
    for k in range(num_iters):
        Q, R = QR_decomposition(A)
        A = R @ Q

    return np.diag(A)


def main():
    A = np.random.rand(3, 3)
    A = A + A.T

    print(np.sort(QR(A)))
    print(np.sort(np.linalg.eigvals(A)))


main()