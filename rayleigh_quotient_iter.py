import numpy as np


def generate_hermite_matrix(n):
    """Generate the n x n real symmetric matrix A used in the example."""
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = 2 * i
            elif i == j + 1 or i == j - 1:
                A[i, j] = -1
    return A


def rayleigh(A, tol, shift, x0):
    """Approximate an eigenvector from initial vector x0 and spectral shift."""
    x = x0 / np.linalg.norm(x0)
    y = np.linalg.solve((A - shift * np.eye(A.shape[0])), x)
    projection = y.T.dot(x)
    shift = shift + 1 / projection
    relative_residual = np.linalg.norm(y - projection * x) / np.linalg.norm(y)

    while relative_residual > tol:
        x = y / np.linalg.norm(y)
        y = np.linalg.solve((A - shift * np.eye(A.shape[0])), x)
        projection = y.T.dot(x)
        shift = shift + 1 / projection
        relative_residual = np.linalg.norm(y - projection * x) / np.linalg.norm(y)

    return x


def main():
    A = generate_hermite_matrix(10)
    x0 = np.random.rand(A.shape[1])
    shift = 1
    tol = 1e-6
    x = rayleigh(A, tol, shift, x0)
    eigval = np.dot(np.dot(x.T, A), x)

    print(eigval)
    print(np.linalg.eig(A)[0])


if __name__ == "__main__":
    main()
