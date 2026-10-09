import numpy as np

from _validation import iteration_limit, matrix, spectral_shift, tolerance, vector
from lanczos import generate_hermite_matrix


def rayleigh(A, tol, shift, x0, max_iter=1000):
    """Return a converged eigenvector by Rayleigh quotient iteration.

    A must be symmetric or Hermitian. Use the supplied shift for the first
    solve, then x* A x for subsequent shifts. Accept only a small eigenpair
    residual relative to ||A||. RQI is not globally convergent; exhaustion
    raises RuntimeError.
    """
    A = matrix(A, symmetric=True)
    tol, max_iter = tolerance(tol), iteration_limit(max_iter)
    shift = spectral_shift(shift)
    x = vector(x0, len(A))
    x /= np.linalg.norm(x)
    scale = max(np.linalg.norm(A, ord=np.inf), np.finfo(float).tiny)
    for k in range(max_iter + 1):
        eigval = np.vdot(x, A @ x).real
        if np.linalg.norm(A @ x - eigval * x) <= tol * scale:
            return x
        if k == max_iter:
            break
        shifted = A - shift * np.eye(len(A))
        try:
            y = np.linalg.solve(shifted, x)
        except np.linalg.LinAlgError as exc:
            # An exact eigenvalue shift has a nullspace: recover its vector.
            _, _, Vh = np.linalg.svd(shifted)
            candidate = Vh[-1].conj()
            if np.linalg.norm(A @ candidate - shift * candidate) <= tol * scale:
                return candidate
            raise RuntimeError("Rayleigh iteration encountered a singular solve") from exc
        y_norm = np.linalg.norm(y)
        if not np.isfinite(y_norm) or y_norm == 0:
            raise FloatingPointError("Rayleigh iteration produced an invalid vector")
        x = y / y_norm
        shift = np.vdot(x, A @ x).real
    raise RuntimeError("Rayleigh iteration did not converge within max_iter")


def main():
    A = generate_hermite_matrix(5)
    x0 = np.array([1., .2, .1, .05, .01])
    x = rayleigh(A, tol=1e-10, shift=-.5, x0=x0)
    eigval = np.vdot(x, A @ x).real
    print(eigval)
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
