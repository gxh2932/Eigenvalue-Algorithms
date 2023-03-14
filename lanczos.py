import numpy as np
from bisection import sturm_bisection, gershgorin_bound


# Function to tri-diagonalize a matrix
def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


def lanczos(A):
    v0 = np.zeros(A.shape[1])
    v0.fill(1.)
    v1 = v0 / np.linalg.norm(v0)

    # First iteration steps
    x, y = [], []
    n = A.shape[1]
    v2, beta = 0.0, 0.0

    for i in range(n):
        # Iteration steps
        w_prime = np.dot(A, v1)
        conj = np.matrix.conjugate(w_prime)
        alpha = np.dot(conj, v1)
        w = w_prime - alpha * v1 - beta * v2
        beta = np.linalg.norm(w)
        x.append(np.linalg.norm(alpha))

        # Reset
        if i < (n-1):
            y.append(beta)
        v2 = v1
        v1 = w/beta

    return tridiag(y, x, y)


def generate_hermite_matrix(n):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = 2 * i
            elif i == j + 1 or i == j - 1:
                A[i, j] = -1
    return A


def main():
    A = generate_hermite_matrix(10)
    T = lanczos(A)

    print(sorted(np.linalg.eig(T)[0]))
    print(sorted(np.linalg.eig(A)[0]))


main()