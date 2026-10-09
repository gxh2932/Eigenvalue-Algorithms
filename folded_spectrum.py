import numpy as np

from _validation import iteration_limit, matrix, spectral_shift, tolerance
from inverse_iter import inverse_iteration


def folded_spectrum_iteration(A, shift, max_iter=1000, tol=1e-10):
    """Return an eigenvector near a real shift for symmetric or Hermitian A.

    Invert the regularized folded matrix (A - shift I)^2. A Rayleigh-Ritz
    step in span{x, A x} resolves mixtures of equally close eigenvalues.
    tol controls the original eigenpair residual relative to ||A||.
    """
    A = matrix(A, symmetric=True)
    shift, tol = spectral_shift(shift), tolerance(tol)
    max_iter = iteration_limit(max_iter)
    if np.iscomplexobj(shift) and np.imag(shift) != 0:
        raise ValueError("Folded spectrum requires a real shift")
    shift = np.real(shift)
    n = len(A)
    x = np.random.rand(n)
    x /= np.linalg.norm(x)
    centered = A - shift * np.eye(n)
    fold_scale = np.linalg.norm(centered, ord=np.inf)
    if fold_scale == 0:
        return x
    centered /= fold_scale
    B = centered @ centered
    # A positive regularization permits an exact eigenvalue target.
    regularization = max(1e-12, n * np.finfo(float).eps)
    shifted = B + regularization * np.eye(n)
    scale = max(np.linalg.norm(A, ord=np.inf), np.finfo(float).tiny)
    for k in range(max_iter):
        y = np.linalg.solve(shifted, x)
        y_norm = np.linalg.norm(y)
        if not np.isfinite(y_norm) or y_norm == 0:
            raise FloatingPointError("Folded spectrum produced an invalid vector")
        x = y / y_norm
        # The folded eigenspace can contain both shift - d and shift + d.
        basis = np.column_stack((x, (A @ x) / scale))
        Q, singular_values, _ = np.linalg.svd(basis, full_matrices=False)
        rank = max(1, int(np.sum(singular_values > tol * singular_values[0])))
        Q = Q[:, :rank]
        eigvals, vectors = np.linalg.eigh(Q.conj().T @ A @ Q)
        index = np.argmin(np.abs(eigvals - shift))
        candidate = Q @ vectors[:, index]
        if np.linalg.norm(A @ candidate - eigvals[index] * candidate) <= tol * scale:
            return candidate
    raise RuntimeError("Folded spectrum did not converge within max_iter")


def main():
    A = np.diag([1., 2., 2.9, 4.])
    x = folded_spectrum_iteration(A, shift=3.)
    eigval = np.vdot(x, A @ x).real
    print(eigval)
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
