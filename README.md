# Eigenvalue Algorithms
Most of these methods have Wikipedia pages. Those that don't I've added a brief description for below.

## Homotopy Method

This method aims to find the eigenpairs (eigenvectors and their corresponding eigenvalues) of a matrix $A$ by solving a system of ODEs derived from a homotopy function. The method starts with the initial conditions based on the standard unit vectors and the diagonal entries of a random diagonal matrix $D$, and gradually transforms the problem into the eigenpairs of matrix $A$ as the homotopy parameter $t$ varies from $0$ to $1$. At the end of the process, the solutions represent the eigenpairs for matrix $A$.

The homotopy-based ODE system can be represented as:

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

with initial conditions $x(0) = e_i$ and $\mu(0) = d_i$ for $i=1,\ldots,n$. As you may infer from the definition of the $i$ variable, the ODE system is solved $n$ times with the solution of each iteration corresponding to the $i\textrm{th}$ eigenpair. Note $x$ represents eigenvector and $\mu$ represents eigenvalue.
