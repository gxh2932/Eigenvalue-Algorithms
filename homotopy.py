import numpy as np
from scipy.integrate import solve_ivp

from _validation import iteration_limit, matrix, tolerance

# Reference: https://www.sciencedirect.com/science/article/pii/0024379588900158


def f(t, y, A, D):
    """Return dy/dt for a real symmetric eigenpair path y = [x, eigval]."""
    A = matrix(A, symmetric=True, real=True)
    D = matrix(D, symmetric=True, real=True)
    y = np.asarray(y)
    if D.shape != A.shape or y.shape != (len(A) + 1,):
        raise ValueError("D and the ODE state must match the dimension of A")
    if np.iscomplexobj(y) or not np.all(np.isfinite(y)) or not np.isfinite(t):
        raise ValueError("The ODE state and t must be finite and real")
    x = y[:-1, None]
    eigval = y[-1]
    if np.linalg.norm(x) == 0:
        raise ValueError("The eigenvector state must be nonzero")
    system_matrix = np.block(
        [[eigval * np.eye(len(A)) - (D + t * (A - D)), x], [x.T, np.zeros((1, 1))]]
    )
    rhs = np.concatenate(((A - D) @ x[:, 0], [0.]))
    try:
        dy_dt = np.linalg.solve(system_matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError(f"Homotopy path is singular at t={t:g}") from exc
    if not np.all(np.isfinite(dy_dt)):
        raise FloatingPointError("Homotopy derivative is not finite")
    return dy_dt


def homotopy_eigenpairs(A, D=None, tol=1e-8, max_iter=10000):
    """Return eigenvalues and eigenvector columns for a real symmetric path.

    D must be diagonal with distinct entries. Each ODE path is limited to
    max_iter right-hand-side evaluations. Crossings or failed integration
    raise RuntimeError; endpoint residuals and orthogonality are verified.
    General nonsymmetric matrices require a different continuation method.
    """
    A = matrix(A, symmetric=True, real=True)
    tol, max_iter = tolerance(tol), iteration_limit(max_iter)
    n = len(A)
    scale = max(np.linalg.norm(A, ord=np.inf), np.finfo(float).tiny)
    if D is None:
        D = np.diag(np.linspace(-scale, scale, n))
    D = matrix(D, real=True)
    if D.shape != A.shape or np.any(D != np.diag(np.diag(D))):
        raise ValueError("D must be a diagonal matrix matching A")
    d = np.diag(D)
    if len(np.unique(d)) != n:
        raise ValueError("D must have distinct diagonal entries")
    # For diagonal A the eigenpairs are known, including repeated entries.
    if np.all(A == np.diag(np.diag(A))):
        order = np.argsort(np.diag(A))
        return np.diag(A)[order].copy(), np.eye(n)[:, order]
    eigvals, Q = np.empty(n), np.empty((n, n))
    for i in range(n):
        y0 = np.concatenate((np.eye(n)[:, i], [d[i]]))
        n_eval = 0

        def rhs(t, y):
            nonlocal n_eval
            n_eval += 1
            if n_eval > max_iter:
                raise RuntimeError(f"Homotopy path {i} exceeded max_iter evaluations")
            return f(t, y, A, D)

        solution = solve_ivp(rhs, (0., 1.), y0, t_eval=[1.],
                             rtol=tol * .01, atol=tol * .01)
        if not solution.success or np.size(solution.y) == 0:
            raise RuntimeError(f"Homotopy path {i} failed: {solution.message}")
        x = solution.y[:-1, -1]
        x /= np.linalg.norm(x)
        eigval = float(x @ A @ x)
        if not np.all(np.isfinite(x)) or np.linalg.norm(A @ x - eigval * x) > tol * scale:
            raise RuntimeError(f"Homotopy path {i} failed its eigenpair residual check")
        eigvals[i], Q[:, i] = eigval, x
    if np.linalg.norm(Q.T @ Q - np.eye(n), ord=np.inf) > 50 * tol:
        raise RuntimeError("Homotopy paths did not produce an orthonormal eigenbasis")
    order = np.argsort(eigvals)
    return eigvals[order], Q[:, order]


def main():
    A = np.array([[2., .2, .1], [.2, 4., .3], [.1, .3, 6.]])
    eigvals, Q = homotopy_eigenpairs(A)
    print(eigvals)
    print(np.linalg.eigvalsh(A))


if __name__ == "__main__":
    main()
