import numpy as np

from _validation import iteration_limit, matrix, tolerance, vector


def conjugate_gradient(A, b, x0, tol=1e-6, max_iter=1000):
    """Solve real symmetric positive-definite A x = b; tol is absolute.

    Return an already converged initial guess immediately. Report detected
    loss of positive definiteness or exhaustion instead of returning NaNs.
    """
    A = matrix(A, symmetric=True, real=True)
    b = vector(b, len(A), nonzero=False)
    x = vector(x0, len(A), nonzero=False)
    if np.iscomplexobj(b) or np.iscomplexobj(x):
        raise ValueError("This routine requires real vectors")
    tol, max_iter = tolerance(tol), iteration_limit(max_iter)
    r = b - A @ x
    if np.linalg.norm(r) <= tol:
        return x
    p = r.copy()
    r_squared = r @ r
    for k in range(max_iter):
        A_p = A @ p
        curvature = p @ A_p
        if not np.isfinite(curvature) or curvature <= 0:
            raise ValueError("Conjugate gradient requires positive-definite A")
        alpha = r_squared / curvature
        x += alpha * p
        r -= alpha * A_p
        r_squared_next = r @ r
        if np.sqrt(r_squared_next) <= tol:
            # Check the true residual before accepting recursive convergence.
            r = b - A @ x
            if np.linalg.norm(r) <= tol:
                return x
            p = r.copy()
            r_squared = r @ r
            continue
        beta = r_squared_next / r_squared
        p = r + beta * p
        r_squared = r_squared_next
    raise RuntimeError("Conjugate gradient did not converge within max_iter")


def main():
    A = np.array([[4., 1., 0.], [1., 3., 1.], [0., 1., 2.]])
    b = np.array([1., 2., 3.])
    x0 = np.zeros(3)
    print(conjugate_gradient(A, b, x0))
    print(np.linalg.solve(A, b))


if __name__ == "__main__":
    main()
