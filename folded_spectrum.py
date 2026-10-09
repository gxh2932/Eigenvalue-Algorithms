# used Bing bot to find this


import numpy as np


def inverse_iteration(A, max_iter: int, shift: float):
    """Apply inverse iteration to the matrix passed by the folded example."""
    x = np.random.rand(A.shape[1])

    for k in range(max_iter):
        y = np.linalg.solve(A - shift * np.eye(A.shape[0])**2, x)

        y_norm = np.linalg.norm(y)

        x = y / y_norm

    return x


def main():
    A = np.random.rand(3, 3)

    shift = 3

    # Transform eigenvalues using the spectral shift
    B = (shift * np.eye(A.shape[0]) - A) @ (shift * np.eye(A.shape[0]) + A)

    # note that B has same eigenvectors as A but different eigenvalues
    x = inverse_iteration(B, 100, shift)
    eigval = np.dot(np.dot(x.T, A), x)

    print(eigval)
    print(np.linalg.eig(A)[0])


if __name__ == "__main__":
    main()
