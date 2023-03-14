# Eigenvalue-Algorithms

## Homotopy Method for Finding Eigenpairs

This Python code aims to find the eigenpairs (eigenvalues and their corresponding eigenvectors) of a matrix A using the homotopy method by solving a system of ODEs. The method starts with the initial conditions based on the standard unit vectors and the eigenvalues of a random diagonal matrix D, and gradually transforms the problem into the eigenpairs of matrix A as the homotopy parameter `t` varies from 0 to 1. At the end of the process, the solutions represent the eigenvalues and their corresponding eigenvectors for matrix A.

### Code Explanation

1. The function `f(t, y, A, D)` represents the right-hand side of the homotopy-based ODE system, with `y` being a vector containing both the variables $x$ (eigenvector) and $\mu$ (eigenvalue).
2. The homotopy-based ODE system is designed such that the solutions of the system at `t=1` will correspond to the eigenpairs of matrix A. The system is solved for $\frac{dy}{dt}$ using `numpy.linalg.solve(M, b)`.
3. The `main` function defines the matrices A and D, with A being the target matrix whose eigenpairs need to be found, and D being a random diagonal matrix.
4. The initial conditions for the ODE system are defined as $x(0) = e_i$, which is the $i$th standard unit vector, and $\mu(0)$ is set to the $i$th diagonal entry of matrix D (the $i$th eigenvalue of D).
5. The ODE system is solved for each initial condition using `scipy.integrate.solve_ivp`. The time span for solving the ODEs is defined as `(0, 1)`. The homotopy parameter `t` varies from 0 to 1, representing the transformation from matrix D to matrix A. The solutions at `t=1` correspond to the eigenpairs of matrix A.
6. The solutions for each $i$ are printed, representing the eigenpairs (eigenvalue `sol.y[-1, -1]` and the corresponding eigenvector `sol.y[:-1, -1]`) of matrix A.
