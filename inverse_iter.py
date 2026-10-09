import numpy as np

from _validation import iteration_limit, matrix, spectral_shift


def inverse_iteration(A, max_iter: int, shift: float):
    """Return a normalized vector after max_iter shifted inverse steps.

    A - shift I must be nonsingular. Convergence requires an isolated
    closest eigenvalue and a starting component in its eigenspace.
    """
    A = matrix(A)
    max_iter, shift = iteration_limit(max_iter), spectral_shift(shift)
    x = np.random.rand(A.shape[1])
    x /= np.linalg.norm(x)
    shifted = A - shift * np.eye(len(A))
    for k in range(max_iter):
        try:
            y = np.linalg.solve(shifted, x)
        except np.linalg.LinAlgError as exc:
            raise np.linalg.LinAlgError(
                "Inverse iteration requires a nonsingular A - shift I"
            ) from exc
        y_norm = np.linalg.norm(y)
        if not np.isfinite(y_norm) or y_norm == 0:
            raise FloatingPointError("Inverse iteration produced an invalid vector")
        x = y / y_norm
    return x


def main():
    A = np.diag([-3., 1., 4.])
    x = inverse_iteration(A, 100, shift=1.1)
    eigval = np.vdot(x, A @ x)
    print(eigval)
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
