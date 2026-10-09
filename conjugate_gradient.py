import numpy as np


def conjugate_gradient(A, b, x0, tol=1e-6, max_iter=1000):
    """Approximate x in A x = b, starting from x0."""
    # Initialize variables
    x = x0
    r = b - A @ x
    p = r
    r_norm = np.linalg.norm(r)

    # Iterate until convergence or maximum iterations
    for k in range(max_iter):
        A_p = A @ p
        alpha = r_norm ** 2 / (p @ A_p)
        x = x + alpha * p
        r = r - alpha * A_p
        r_norm_next = np.linalg.norm(r)
        if r_norm_next < tol:
            break
        beta = r_norm_next ** 2 / r_norm ** 2
        p = r + beta * p
        r_norm = r_norm_next

    return x


def main():
    A = np.array([[1, 2, 3], [2, 5, 6], [3, 6, 9]])
    b = np.array([1, 2, 3])
    x0 = np.array([1, 1, 1])
    x = conjugate_gradient(A, b, x0)
    print(x)
    print(np.linalg.solve(A, b))


if __name__ == "__main__":
    main()
