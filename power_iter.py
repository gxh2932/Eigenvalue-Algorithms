import numpy as np

from _validation import iteration_limit, matrix


def power_iteration(A, max_iter: int):
    """Return a normalized vector after max_iter power steps.

    Convergence requires an isolated dominant eigenvalue in magnitude and a
    nonzero initial component in its eigenspace. This fixed-count routine
    returns an approximation; it does not certify convergence.
    """
    A = matrix(A)
    max_iter = iteration_limit(max_iter)
    x = np.random.rand(A.shape[1])
    x /= np.linalg.norm(x)
    for k in range(max_iter):
        y = A @ x
        y_norm = np.linalg.norm(y)
        if y_norm == 0:
            return x  # A x = 0: x is already an eigenvector.
        if not np.isfinite(y_norm):
            raise FloatingPointError("Power iteration exceeded floating-point range")
        x = y / y_norm
    return x


def main():
    A = np.diag([-5., 2., 1.])
    x = power_iteration(A, 100)
    eigval = np.vdot(x, A @ x)
    print(eigval)
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
