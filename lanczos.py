import numpy as np
from bisection import sturm_bisection, gershgorin_bound


def tridiag(e_lower, d, e_upper, lower_offset=-1, diagonal_offset=0, upper_offset=1):
    """Construct T from its lower, main, and upper diagonal entries."""
    return (
        np.diag(e_lower, lower_offset)
        + np.diag(d, diagonal_offset)
        + np.diag(e_upper, upper_offset)
    )


def lanczos(A):
    """Return the tridiagonal matrix T from the Lanczos recurrence for A."""
    x0 = np.zeros(A.shape[1])
    x0.fill(1.)
    q = x0 / np.linalg.norm(x0)

    # First iteration steps
    d, e = [], []
    n = A.shape[1]
    q_prev, beta = 0.0, 0.0

    for k in range(n):
        # Iteration steps
        y = np.dot(A, q)
        y_conj = np.matrix.conjugate(y)
        alpha = np.dot(y_conj, q)
        r = y - alpha * q - beta * q_prev
        beta = np.linalg.norm(r)
        d.append(np.linalg.norm(alpha))

        # Reset
        if k < (n - 1):
            e.append(beta)
        q_prev = q
        q = r / beta

    return tridiag(e, d, e)


def generate_hermite_matrix(n):
    """Generate the n x n real symmetric matrix A used in the example."""
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


if __name__ == "__main__":
    main()
