# Eigenvalue Algorithms
Educational implementations of eigenvalue algorithms, plus a conjugate-gradient linear solver. The routines use the notation below and check their input domains and numerical stopping conditions.

## Running and testing

Python 3.10 or newer is required.

```sh
python -m pip install -r requirements.txt
python power_iter.py
python -m unittest discover -s tests -v
```

The tests compare eigenvalues, eigenpair residuals, orthogonality, and linear-system solutions against NumPy. They include the original failing cases, signed and repeated spectra, zero pivots, complex inputs where supported, input preservation, iteration exhaustion, and every executable demonstration. GitHub Actions runs the suite on Python 3.10 and 3.12.

## Input domains and stopping conditions

| Routine | Supported input / behavior |
| --- | --- |
| QR eigenvalue iteration and Jacobi | Real symmetric matrices; return eigenvalues only after relative convergence. QR uses NumPy's Householder QR decomposition, Wilkinson shifts, and deflation. |
| Arnoldi | Real or complex square matrices and a nonzero starting vector. Return the active basis and Hessenberg matrix on Krylov breakdown. |
| Lanczos | Symmetric or Hermitian matrices; signed recurrence coefficients and reorthogonalization. On breakdown, the smaller returned matrix represents the reached invariant subspace, including only its eigenvalue multiplicities. |
| Bisection / divide and conquer | Real symmetric tridiagonal matrices, including scalar matrices and zero off-diagonal entries. Divide and conquer preserves its input and uses NumPy's symmetric eigensolver for the rank-one merge. |
| Conjugate gradient | Real symmetric positive-definite systems. Return an already solved initial guess immediately; reject detected nonpositive curvature. |
| Power / inverse iteration | Real or complex square matrices; fixed iteration counts return approximate vectors. Power convergence requires an isolated dominant eigenvalue in magnitude. Inverse convergence requires an isolated closest eigenvalue and a nonsingular shifted matrix. Both need an initial component in the desired eigenspace. |
| Rayleigh quotient iteration | Symmetric or Hermitian matrices. Check the original eigenpair residual and update shifts from the Rayleigh quotient. Convergence is local; a cycling start reports exhaustion. |
| Folded spectrum | Symmetric or Hermitian matrices and a real target. Use $(A-\sigma I)^2$, positive regularization for exact targets, and a small Rayleigh-Ritz step to resolve equally close eigenvalues. |
| Laguerre | Finite nonconstant polynomials, supplied as `np.poly1d` or coefficient arrays in descending order. Real and complex roots are supported. Forming characteristic-polynomial coefficients and repeatedly deflating high-degree or clustered polynomials can be ill-conditioned; use tridiagonal bisection for those eigenvalue problems. |
| Homotopy | Real symmetric matrices and a diagonal starting matrix with distinct entries. Check integration success, endpoint residuals, and basis orthogonality. Singular paths or crossings can prevent continuation and report an error. General nonsymmetric or complex targets require a different continuation method. |

`tol` is an absolute interval width in bisection and an absolute residual norm in conjugate gradient. QR, Jacobi, Rayleigh, folded spectrum, and homotopy use residual or off-diagonal tolerances relative to the matrix norm. Arnoldi and Lanczos use relative breakdown tolerances. Laguerre uses a coefficient-scaled polynomial residual. Bounded convergence routines raise `RuntimeError` on exhaustion; invalid inputs raise `ValueError`, and singular inverse solves or arithmetic failures are reported explicitly.

## Notation

The code and descriptions use the following conventions for the eigenvalue problem $Ax = \lambda x$.

| Quantity | Mathematical notation | Python name |
| --- | --- | --- |
| Input matrix | $A$ | `A` |
| Symmetric tridiagonal matrix | $T$ | `T` |
| Diagonal matrix | $D$ | `D` |
| Orthogonal basis or eigenvector matrix | $Q$ | `Q` (vectors are columns) |
| Upper triangular QR factor | $R$ | `R` |
| Upper Hessenberg Arnoldi matrix | $H$ | `H` |
| Eigenvector or current vector iterate | $x$ | `x` |
| Initial vector | $x_0$ | `x0` |
| Unnormalized next vector | $y$ | `y` |
| Eigenvalue / collection of eigenvalues | $\lambda$ / $\lambda_i$ | `eigval` / `eigvals` |
| Spectral shift | $\sigma$ | `shift` |
| Matrix dimension | $n$ | `n` |
| Number of Krylov steps | $m$ | `m` |
| Iteration index / matrix indices | $k$ / $i, j$ | `k` / `i`, `j` |
| Convergence tolerance | $\mathrm{tol}$ | `tol` |
| Iteration limit or fixed iteration count | $k_{\max}$ | `max_iter` |
| Completed iteration count | — | `n_iter` |
| Diagonal / off-diagonal entries of $T$ | $d_i$ / $e_i$ | `d` / `e` |
| Lower / upper eigenvalue bounds | $\ell$ / $u$ | `lower` / `upper` |
| Polynomial variable / initial root guess | $z$ / $z_0$ | `z` / `z0` |

Python arrays use zero-based indices. Mathematical subscripts and the `index` argument to `sturm_bisection` use one-based indices. Rectangular QR decomposition uses `m` rows and `n` columns. Algorithm-specific scalars such as `alpha` and `beta` retain their usual local meanings in conjugate gradient and Lanczos; a Lanczos basis vector is `q`, its predecessor is `q_prev`, and a residual is `r`. In homotopy, `y` is the combined ODE state `[x, eigval]`.

### Calling the routines

Existing function and file names are retained. Keyword arguments now follow the conventions above; the original positional argument order is retained, with new optional controls appended. Arnoldi now returns only the active columns on breakdown; divide and conquer always returns a one-dimensional eigenvalue array and preserves the input matrix.

| Routine | Current signature | Renamed arguments |
| --- | --- | --- |
| QR iteration | `QR(A, max_iter=1000, tol=1e-12)` | `num_iters` → `max_iter`; adds relative `tol` |
| Power iteration | `power_iteration(A, max_iter)` | `num_iterations` → `max_iter` |
| Inverse iteration (both modules) | `inverse_iteration(A, max_iter, shift)` | `num_iterations` → `max_iter`, `mu` → `shift` |
| Arnoldi iteration | `arnoldi_iteration(A, x0, m, tol=1e-12)` | `b` → `x0`, `n` → `m`; `tol` exposes the previous fixed tolerance |
| Lanczos | `lanczos(A, m=None, tol=1e-12, x0=None)` | Adds optional Krylov dimension, breakdown tolerance, and starting vector |
| Jacobi iteration | `jacobi_eigenvalue_algorithm(A, tol=1e-10, max_iter=10000)` | `tolerance` → `tol`; adds iteration limit |
| Rayleigh iteration | `rayleigh(A, tol, shift, x0, max_iter=1000)` | `epsilon` → `tol`, `mu` → `shift`, `x` → `x0`; adds iteration limit |
| Sturm sequence | `sturm_evaluate(z, d, e)` | `t` → `z`, `a` → `d`, `b` → `e` |
| Sturm bisection | `sturm_bisection(index, d, e, lower, upper, tol=1e-6, max_iter=1000)` | `k` → `index`, `a` → `d`, `b` → `e`, `alpha` → `lower`, `beta` → `upper`; adds tolerance and iteration controls |
| Gershgorin bounds (both modules) | `gershgorin_bound(d, e)` | `a` → `d`, `b` → `e` |
| Characteristic polynomial | `construct_characteristic_polynomial(d, e)` | `a` → `d`, `b` → `e` |
| Laguerre iteration | `laguerre_method(p, z0, tol=1e-6, max_iter=100)` | `x0` → `z0`, `epsilon` → `tol` |
| Divide and conquer | `div_conq(T)`, `partition(T)` | `A` → `T` |
| Folded spectrum | `folded_spectrum_iteration(A, shift, max_iter=1000, tol=1e-10)` | New callable routine for the corrected folded transformation |
| Homotopy | `homotopy_eigenpairs(A, D=None, tol=1e-8, max_iter=10000)` | New callable routine; `max_iter` caps RHS evaluations per path |
| Tridiagonal construction | `tridiag(e_lower, d, e_upper, lower_offset=-1, diagonal_offset=0, upper_offset=1)` | `a`, `b`, `c` → `e_lower`, `d`, `e_upper`; `k1`, `k2`, `k3` → the corresponding offsets |

For example, use `inverse_iteration(A, max_iter=100, shift=1.0)` and `laguerre_method(p, z0=0.0, tol=1e-6)`. Existing calls that use the old keyword names must be updated. `conjugate_gradient(A, b, x0, tol=1e-6, max_iter=1000)` keeps `b` for the right-hand side of $Ax=b$.

The algorithms require NumPy; homotopy also requires SciPy. Run a demonstration with, for example, `python power_iter.py`. Importing a module makes its routines available without running its demonstration.

## Bisection Method
Essentially a variation of the classical bisection method, which is a root-finding algorithm that works by repeatedly dividing an interval in half and then selecting the subinterval where the eigenvalue is guaranteed to lie. This algorithm finds the roots of the characteristic polynomial for a $n \times n$ symmetric tridiagonal matrix, using the Gershgorin Circle Theorem (a personal favorite of mine) and Sturm sequences.

The method works by first computing a lower bound $\ell$ and upper bound $u$ for the eigenvalues of $T$, using the Gershgorin Circle Theorem. Once a bound is found, we iterate through each eigenvalue index $i$ from $1$ to $n$. For each eigenvalue index $i$, we apply the bisection method augmented with Sturm sequences to find the $i$\-th smallest eigenvalue $\lambda_i$ of $T$.

The Sturm routine counts eigenvalues strictly greater than $z$ using the pivots of an $LDL^T$ factorization of $T-zI$. Bisection uses `n - count`, the number of eigenvalues less than or equal to $z$. Small negative replacements for zero pivots handle equality consistently, and zero off-diagonal entries start independent blocks. This avoids multiplying characteristic-polynomial values, which can overflow. The equivalent determinant sequences $q_i(z)=\det(zI-T_i)$ and $p_i(z)=\det(T_i-zI)$ satisfy $q_i(z)=(-1)^i p_i(z)$, where $T_i$ is the leading principal block. Laguerre uses the latter convention.

For more info refer to _Numerical Methods for Eigenvalue Problems_ (2012) by Steffen Börm.


## Homotopy Method

This implementation finds eigenpairs of a real symmetric $n \times n$ matrix $A$ by solving a system of ODEs derived from a homotopy function. It starts with a standard unit vector and an entry of a diagonal matrix $D$ with distinct entries, and gradually transforms the problem into an eigenpair of $A$ as the parameter $t$ changes from $0$ to $1$. The implementation verifies each endpoint and reports paths that fail. The homotopy-based ODE system is:

```math
\begin{bmatrix}
    \lambda I - [D + t(A - D)] & x \\
    x^* & 0
\end{bmatrix}
\begin{bmatrix}
    \dot{x} \\
    \dot{\lambda}
\end{bmatrix}
=
\begin{bmatrix}
    (A - D)x \\
    0
\end{bmatrix},
```

with initial conditions $x(0) = \mathbf{e}_i$ and $\lambda(0) = d_i$ for $i=1,\ldots,n$, where $\mathbf{e}_i$ is the $i$\-th standard unit vector, $d_i=D_{ii}$, $x$ represents the eigenvector, and $\lambda$ represents the eigenvalue. The superscript $*$ denotes conjugate transpose. The ODE system is solved $n$ times, with each solution corresponding to the $i$\-th eigenpair. In code, the standard unit vector is `e_i` and the eigenvalue is `eigval`.

More info (such as the derivation) can be found here: https://www.sciencedirect.com/science/article/pii/0024379588900158


## Laguerre Iteration

Method for finding the eigenvalues of a $n \times n$ symmetric tridiagonal matrix. The algorithm begins by constructing the characteristic polynomial of the given matrix using the following recurrence relation for symmetric tridiagonal matrices:

```math
p_i(z) = (d_i - z)p_{i-1}(z) - e_{i-1}^2 p_{i-2}(z),
```

where $d_i$ is the $i$\-th diagonal entry and $e_i$ is the entry connecting rows $i$ and $i+1$ of $T$. The initial polynomials are $p_0(z)=1$ and $p_1(z)=d_1-z$. The method then utilizes Laguerre's method, an iterative root-finding technique, to approximate the eigenvalues by finding the roots of $p_n(z)=\det(T-zI)$. Gershgorin bounds provide the initial root guess $z_0=(\ell+u)/2$. Finally, the code iteratively updates the characteristic polynomial by dividing out the factors corresponding to the found eigenvalues.
