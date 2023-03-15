# Eigenvalue Algorithms
Most of these methods have Wikipedia pages. Those that don't I've added a brief description for below.

## Bisection Method
Essentially a variation of the classical bisection method, which is a root-finding algorithm that works by repeatedly dividing an interval in half and then selecting the subinterval where the eigenvalue is guaranteed to lie. This algorithm finds the roots of the characteristic polynomial for a symmetric tridiagonal matrix, using the Gershgorin Circle Theorem (a personal favorite of mine) and Sturm sequences.

The method works by first computing a lower bound alpha and upper bound beta for the eigenvalues of the matrix, using the Gershgorin Circle Theorem. Once a bound is found, we iterate through each eigenvalue index k from 1 to n, where n is the size of the matrix. For each eigenvalue index k, we apply the bisection method augmented with Sturm sequences to find the k-th smallest eigenvalue of the matrix.

When applying the bisection method to find an eigenvalue, we start with the interval found using Gershgorin (you could do this more efficiently by storing the results of the Sturm sequence calculations from prior bisection applications but we don't really care about efficiency here). The midpoint of this interval, gamma, is calculated, and we need to determine whether the eigenvalue we're looking for lies in the subinterval [alpha, gamma] or [gamma, beta]. This decision is based on the number of eigenvalues less than or equal to gamma. Sturm sequences provide an efficient way to compute this information without explicitly finding all the eigenvalues of the matrix. By evaluating the Sturm sequence for the value gamma, we can count the number of eigenvalues less than or equal to gamma. Based on this count, we can decide which subinterval to choose for the next iteration of the bisection method. This process is repeated until convergence, as seen in the classical bisection method.

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
