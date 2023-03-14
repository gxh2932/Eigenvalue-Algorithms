import numpy as np


# real symmetric matrix


def jacobi_eigenvalue_algorithm(A, tolerance=1e-10):
    n = A.shape[0]
    eigenvectors = np.eye(n)
    iterations = 0
    while True:
        # Find maximum off-diagonal element
        max_offdiag = 0
        max_i, max_j = 0, 0
        for i in range(n):
            for j in range(i+1, n):
                if abs(A[i, j]) > max_offdiag:
                    max_offdiag = abs(A[i, j])
                    max_i, max_j = i, j

        if max_offdiag < tolerance:
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
        eigenvectors = np.dot(eigenvectors, J)
        iterations += 1

    # Extract eigenvalues and eigenvectors
    eigenvalues = np.diag(A)

    return eigenvalues, eigenvectors, iterations


def main():
    A = np.random.randn(3,3)
    A = A + A.T

    eigenvalues, eigenvectors, iterations = jacobi_eigenvalue_algorithm(A)

    print(eigenvalues)
    print(np.linalg.eig(A)[0])


main()