# Eigenvalue Algorithms
Most of these methods have Wikipedia pages. Those that don't I've added a brief description for below.

## Bisection Method
Essentially a variation of the classical bisection method, which is a root-finding algorithm that works by repeatedly dividing an interval in half and then selecting the subinterval where the eigenvalue is guaranteed to lie. This algorithm finds the roots of the characteristic polynomial for a $n \times n$ symmetric tridiagonal matrix, using the Gershgorin Circle Theorem (a personal favorite of mine) and Sturm sequences.

The method works by first computing a lower bound alpha and upper bound beta for the eigenvalues of the matrix, using the Gershgorin Circle Theorem. Once a bound is found, we iterate through each eigenvalue index $i$ from $1$ to $n$. For each eigenvalue index $i$, we apply the bisection method augmented with Sturm sequences to find the $i$\-th smallest eigenvalue of the matrix.

For more info refer to _Numerical Methods for Eigenvalue Problems_ (2012) by Steffen Börm.


## Homotopy Method

This method aims to find the eigenpairs (eigenvectors and their corresponding eigenvalues) of a $n \times n$ matrix $A$ by solving a system of ODEs derived from a homotopy function. The method starts with the initial conditions based on a standard unit vector and diagonal entry of a random diagonal matrix $D$, and gradually transforms the problem into an eigenpair of matrix $A$ as the homotopy parameter $t$ changes from $0$ to $1$. The homotopy-based ODE system can be represented as:

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

with initial conditions $x(0) = e_i$ and $\mu(0) = d_i$ for $i=1,\ldots,n$ where $x$ represents eigenvector and $\mu$ represents eigenvalue. As you may infer from the definition of the $i$ variable, the ODE system is solved $n$ times with the solution of each iteration corresponding to the $i$\-th eigenpair.

More info (such as the derivation) can be found here: https://www.sciencedirect.com/science/article/pii/0024379588900158


## Laguere Iteration

Method for finding the eigenvalues of a $n \times n$ symmetric tridiagonal matrix. The algorithm begins by constructing the characteristic polynomial of the given matrix using the following recurrence relation for symmetric tridiagonal matrices:

```math
p_i(x) = (a_i - x)p_{i-1}(x) - b_{i-1}^2 p_{i-2}(x),
```

where $a_i$ is the $i\text{th}$ diagonal element and $b_i$ is the $i$\-th subdiagonal element. It then utilizes Laguerre's method, an iterative root-finding technique, to approximate the eigenvalues of the matrix by finding the roots of the characteristic polynomial. Gershgorin Circle Theorem is employed to provide an initial guess for the eigenvalue search and to bound the search interval. Finally, the code iteratively updates the characteristic polynomial by dividing out the factors corresponding to the found eigenvalues, ensuring that all eigenvalues are approximated.
