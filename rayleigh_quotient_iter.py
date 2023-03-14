import numpy as np


def generate_hermite_matrix(n):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = 2 * i
            elif i == j + 1 or i == j - 1:
                A[i, j] = -1
    return A


def rayleigh(A, epsilon, mu, x):
    x = x / np.linalg.norm(x)
    y = np.linalg.solve((A - mu * np.eye(A.shape[0])), x)
    lambda_val = y.T.dot(x)
    mu = mu + 1 / lambda_val
    err = np.linalg.norm(y - lambda_val * x) / np.linalg.norm(y)

    while err > epsilon:
        x = y / np.linalg.norm(y)
        y = np.linalg.solve((A - mu * np.eye(A.shape[0])), x)
        lambda_val = y.T.dot(x)
        mu = mu + 1 / lambda_val
        err = np.linalg.norm(y - lambda_val * x) / np.linalg.norm(y)

    return x


def main():
    A = generate_hermite_matrix(10)
    x = np.random.rand(A.shape[1])
    mu = 1
    epsilon = 1e-6
    eigenvector = rayleigh(A, epsilon, mu, x)
    eigenvalue = np.dot(np.dot(eigenvector.T, A), eigenvector)

    print(eigenvalue)
    print(np.linalg.eig(A)[0])


main()