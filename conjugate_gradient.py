import numpy as np


def conjugate_gradient(A, b, x0, tol=1e-6, max_iter=1000):
    # Initialize variables
    x = x0
    r = b - A @ x
    p = r
    r_norm = np.linalg.norm(r)

    # Iterate until convergence or maximum iterations
    for i in range(max_iter):
        Ap = A @ p
        alpha = r_norm ** 2 / (p @ Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        r_norm_new = np.linalg.norm(r)
        if r_norm_new < tol:
            break
        beta = r_norm_new ** 2 / r_norm ** 2
        p = r + beta * p
        r_norm = r_norm_new

    return x


def main():
    A = np.array([[1, 2, 3], [2, 5, 6], [3, 6, 9]])
    b = np.array([1, 2, 3])
    x0 = np.array([1, 1, 1])
    x = conjugate_gradient(A, b, x0)
    print(x)
    print(np.linalg.solve(A, b))


main()