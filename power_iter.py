import numpy as np


def power_iteration(A, num_iterations: int):
    # Ideally choose a random vector
    # To decrease the chance that our vector
    # Is orthogonal to the eigenvector
    b_k = np.random.rand(A.shape[1])

    for _ in range(num_iterations):
        # calculate the matrix-by-vector product Ab
        b_k1 = np.dot(A, b_k)

        # calculate the norm
        b_k1_norm = np.linalg.norm(b_k1)

        # re normalize the vector
        b_k = b_k1 / b_k1_norm

    return b_k


def main():
    A = np.random.rand(10, 10)

    eigenvector = power_iteration(A, 100)
    max_eigenvalue = np.dot(np.dot(eigenvector.T, A), eigenvector)

    print(max_eigenvalue)


main()