import numpy as np


def jacobi_eigenvalue_algorithm(A, tol=1e-10):
    """Return eigvals, eigenvector columns Q, and n_iter for real symmetric A."""
    n = A.shape[0]
    Q = np.eye(n)
    n_iter = 0
    while True:
        # Find maximum off-diagonal element
        max_offdiag = 0
        max_i, max_j = 0, 0
        for i in range(n):
            for j in range(i+1, n):
                if abs(A[i, j]) > max_offdiag:
                    max_offdiag = abs(A[i, j])
                    max_i, max_j = i, j

        if max_offdiag < tol:
            break

        # Compute the Jacobi rotation matrix
        theta = 0.5 * np.arctan2(2 * A[max_i, max_j], A[max_i, max_i] - A[max_j, max_j])
        c = np.cos(theta)
        s = np.sin(theta)
        J = np.eye(n)
        J[max_i, max_i] = c
        J[max_j, max_j] = c
        J[max_i, max_j] = -s
        J[max_j, max_i] = s

        # Update the matrix and eigenvectors
        A = np.dot(np.dot(J.T, A), J)
        Q = np.dot(Q, J)
        n_iter += 1

    # Extract eigenvalues and eigenvectors
    eigvals = np.diag(A)

    return eigvals, Q, n_iter


def main():
    A = np.random.randn(3,3)
    A = A + A.T

    eigvals, Q, n_iter = jacobi_eigenvalue_algorithm(A)

    print(eigvals)
    print(np.linalg.eig(A)[0])


if __name__ == "__main__":
    main()
