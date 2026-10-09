import numpy as np


def inverse_iteration(A, max_iter: int, shift: float):
    """Approximate an eigenvector near shift using inverse iteration."""
    x = np.random.rand(A.shape[1])

    for k in range(max_iter):
        y = np.linalg.solve(A - shift * np.eye(A.shape[0]), x)

        y_norm = np.linalg.norm(y)

        x = y / y_norm

    return x


def main():
    A = np.random.rand(10, 10)

    x = inverse_iteration(A, 100, 1)
    eigval = np.dot(np.dot(x.T, A), x)

    print(eigval)
    print(np.linalg.eig(A)[0])


if __name__ == "__main__":
    main()
