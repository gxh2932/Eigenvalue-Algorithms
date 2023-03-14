# used Bing bot to find this


import numpy as np


def inverse_iteration(A, num_iterations: int, mu: float):
    b_k = np.random.rand(A.shape[1])

    for _ in range(num_iterations):
        b_k1 = np.linalg.solve(A - mu * np.eye(A.shape[0])**2, b_k)

        b_k1_norm = np.linalg.norm(b_k1)

        b_k = b_k1 / b_k1_norm

    return b_k


def main():
    A = np.random.rand(3, 3)

    pole = 3

    # fold eigenvalues of A around pole
    B = (pole * np.eye(A.shape[0]) - A) @ (pole * np.eye(A.shape[0]) + A)

    # note that B has same eigenvectors as A but different eigenvalues
    eigenvector = inverse_iteration(B, 100, pole)
    eigenvalue = np.dot(np.dot(eigenvector.T, A), eigenvector)

    print(eigenvalue)
    print(np.linalg.eig(A)[0])


main()