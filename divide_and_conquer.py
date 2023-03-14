import numpy as np

# reference: https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf


def partition(A):
    n = A.shape[0]
    m = n // 2

    A1 = A[:m, :m]
    A1[m-1, m-1] = A1[m-1, m-1] - A[m-1, m]

    A2 = A[m:, m:]
    A2[0, 0] = A2[0, 0] - A[m-1, m]

    T = np.zeros((n, n))
    T[:m, :m] = A1
    T[m:, m:] = A2

    U = np.zeros((n, n))
    U[m-1, m-1] = A[m-1, m]
    U[m-1, m] = A[m-1, m]
    U[m, m-1] = A[m-1, m]
    U[m, m] = A[m-1, m]

    return T, U


def div_conq(A):
    n = A.shape[0]

    if A.shape == (1, 1):
        return A, np.eye(1)
    else:
        T, U = partition(A)
        alpha1, Q1 = div_conq(T[:T.shape[0]//2, :T.shape[0]//2])
        alpha2, Q2 = div_conq(T[T.shape[0]//2:, T.shape[0]//2:])

        D = np.zeros((n, n))
        alpha_1_len = alpha1.shape[0]
        for i in range(n):
            if i < alpha_1_len:
                D[i, i] = alpha1[0]
                alpha1 = alpha1[1:]
            else:
                D[i, i] = alpha2[0]
                alpha2 = alpha2[1:]

        V_top = Q1[-1]
        V_bottom = Q2[0]
        V = np.concatenate((V_top, V_bottom), axis=0)
        V = V.reshape((n, 1))

        ro = A[n//2-1, n//2]
        V = ro * (V @ V.T)

        M = D + V

        alpha, Q_prime = np.linalg.eigh(M)  # kind of cheating

        Q = np.zeros((n, n))
        Q[:n//2, :n//2] = Q1
        Q[n//2:, n//2:] = Q2

        Q = Q @ Q_prime

    return alpha, Q


def symmetric_tridiagonal_matrix(n):
    d = np.random.rand(n)
    e = np.random.rand(n-1)
    T = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
    return T


def main():
    T = symmetric_tridiagonal_matrix(5)
    alpha, Q = div_conq(T)
    print(sorted(alpha))
    print(sorted(np.linalg.eigvals(T)))


main()