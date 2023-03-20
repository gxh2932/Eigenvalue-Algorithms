# Eigenvalue Algorithms
Most of these methods have Wikipedia pages. Those that don't I've added a brief description for below.

## Bisection Method
Essentially a variation of the classical bisection method, which is a root-finding algorithm that works by repeatedly dividing an interval in half and then selecting the subinterval where the eigenvalue is guaranteed to lie. This algorithm finds the roots of the characteristic polynomial for a symmetric tridiagonal matrix, using the Gershgorin Circle Theorem (a personal favorite of mine) and Sturm sequences.

The method works by first computing a lower bound alpha and upper bound beta for the eigenvalues of the matrix, using the Gershgorin Circle Theorem. Once a bound is found, we iterate through each eigenvalue index k from 1 to n, where n is the size of the matrix. For each eigenvalue index k, we apply the bisection method augmented with Sturm sequences to find the k-th smallest eigenvalue of the matrix.

For more info refer to Numerical Methods for Eigenvalue Problems (2012) by Steffen Börm.


## Homotopy Method

This method aims to find the eigenpairs (eigenvectors and their corresponding eigenvalues) of a matrix $A$ by solving a system of ODEs derived from a homotopy function. The method starts with the initial conditions based on a standard unit vector and diagonal entry of a random diagonal matrix $D$, and gradually transforms the problem into an eigenpair of matrix $A$ as the homotopy parameter $t$ changes from $0$ to $1$. The homotopy-based ODE system can be represented as:

```math
\begin{bmatrix}
    \mu I - [D + t(A - D)] & x \\
    x^* & 0
\end{bmatrix}
\begin{bmatrix}
    \dot{x} \\
    \dot{\mu}
\end{bmatrix}
=
\begin{bmatrix}
    (A - D)x \\
    0
\end{bmatrix},
```

with initial conditions $x(0) = e_i$ and $\mu(0) = d_i$ for $i=1,\ldots,n$, where $x$ represents eigenvector and $\mu$ represents eigenvalue. As you may infer from the definition of the $i$ variable, the ODE system is solved $n$ times with the solution of each iteration corresponding to the $i\text{th}$ eigenpair.

More info (such as the derivation) can be found here: https://www.sciencedirect.com/science/article/pii/0024379588900158


## Laguere Iteration

Laguerre's method for finding the eigenvalues of a symmetric tridiagonal matrix. The algorithm begins by constructing the characteristic polynomial of the given matrix using the following recurrence relation for symmetric tridiagonal matrices:

```math
p_n(x) = (a_n - x)p_{n-1}(x) - b_{n-1}^2 p_{n-2}(x).
```

It then utilizes Laguerre's method, an iterative root-finding technique, to approximate the eigenvalues of the matrix by finding the roots of the characteristic polynomial. Gershgorin Circle Theorem is employed to provide an initial guess for the eigenvalue search and to bound the search interval. Finally, the code iteratively updates the characteristic polynomial by dividing out the factors corresponding to the found eigenvalues, ensuring that all eigenvalues are approximated.
