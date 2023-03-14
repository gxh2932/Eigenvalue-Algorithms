import numpy as np


def inverse_iteration(A, num_iterations: int, mu: float):
    b_k = np.random.rand(A.shape[1])

    for _ in range(num_iterations):
        b_k1 = np.linalg.solve(A - mu * np.eye(A.shape[0]), b_k)

        b_k1_norm = np.linalg.norm(b_k1)

        b_k = b_k1 / b_k1_norm

    return b_k


def main():
    A = np.random.rand(10, 10)

    eigenvector = inverse_iteration(A, 100, 1)
    eigenvalue = np.dot(np.dot(eigenvector.T, A), eigenvector)

    print(eigenvalue)
    print(np.linalg.eig(A)[0])


main()