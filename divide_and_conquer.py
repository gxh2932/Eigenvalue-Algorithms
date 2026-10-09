import numpy as np

# reference: https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf


def partition(T):
    """Split tridiagonal T into block diagonal B and a rank-one update."""
    n = T.shape[0]
    split = n // 2

    T1 = T[:split, :split]
    T1[split - 1, split - 1] = T1[split - 1, split - 1] - T[split - 1, split]

    T2 = T[split:, split:]
    T2[0, 0] = T2[0, 0] - T[split - 1, split]

    B = np.zeros((n, n))
    B[:split, :split] = T1
    B[split:, split:] = T2

    rank_one_update = np.zeros((n, n))
    rank_one_update[split - 1, split - 1] = T[split - 1, split]
    rank_one_update[split - 1, split] = T[split - 1, split]
    rank_one_update[split, split - 1] = T[split - 1, split]
    rank_one_update[split, split] = T[split - 1, split]

    return B, rank_one_update


def div_conq(T):
    """Compute eigenvalues and eigenvector columns Q of tridiagonal T."""
    n = T.shape[0]

    if T.shape == (1, 1):
        return T, np.eye(1)
    else:
        B, rank_one_update = partition(T)
        eigvals1, Q1 = div_conq(B[:B.shape[0]//2, :B.shape[0]//2])
        eigvals2, Q2 = div_conq(B[B.shape[0]//2:, B.shape[0]//2:])

        D = np.zeros((n, n))
        n1 = eigvals1.shape[0]
        for i in range(n):
            if i < n1:
                D[i, i] = eigvals1[0]
                eigvals1 = eigvals1[1:]
            else:
                D[i, i] = eigvals2[0]
                eigvals2 = eigvals2[1:]

        q_top = Q1[-1]
        q_bottom = Q2[0]
        v = np.concatenate((q_top, q_bottom), axis=0)
        v = v.reshape((n, 1))

        rho = T[n//2 - 1, n//2]
        rank_one_update = rho * (v @ v.T)

        secular_matrix = D + rank_one_update

        eigvals, Q_update = np.linalg.eigh(secular_matrix)  # kind of cheating

        Q = np.zeros((n, n))
        Q[:n//2, :n//2] = Q1
        Q[n//2:, n//2:] = Q2

        Q = Q @ Q_update

    return eigvals, Q


def symmetric_tridiagonal_matrix(n):
    """Generate an n x n real symmetric tridiagonal matrix T."""
    d = np.random.rand(n)
    e = np.random.rand(n-1)
    T = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
    return T


def main():
    T = symmetric_tridiagonal_matrix(5)
    eigvals, Q = div_conq(T)
    print(sorted(eigvals))
    print(sorted(np.linalg.eigvals(T)))


if __name__ == "__main__":
    main()
