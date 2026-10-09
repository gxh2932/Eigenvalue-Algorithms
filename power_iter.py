import numpy as np


def power_iteration(A, max_iter: int):
    """Approximate a dominant eigenvector of A using power iteration."""
    # Ideally choose a random vector
    # To decrease the chance that our vector
    # Is orthogonal to the eigenvector
    x = np.random.rand(A.shape[1])

    for k in range(max_iter):
        # Calculate the matrix-vector product A x
        y = np.dot(A, x)

        # calculate the norm
        y_norm = np.linalg.norm(y)

        # re normalize the vector
        x = y / y_norm

    return x


def main():
    A = np.random.rand(10, 10)

    x = power_iteration(A, 100)
    eigval = np.dot(np.dot(x.T, A), x)

    print(eigval)


if __name__ == "__main__":
    main()
