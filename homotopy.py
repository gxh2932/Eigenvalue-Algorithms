import numpy as np
import scipy.integrate as spi

# reference: https://www.sciencedirect.com/science/article/pii/0024379588900158


# Define the right-hand side of the ODE system
def f(t, y, A, D):
    """Return dy/dt for the eigenpair state y = [x, eigval]."""
    y = y.reshape(y.shape[0], 1)

    x = y[:-1]
    eigval = y[-1]

    # Compute the augmented system matrix and right-hand side
    system_matrix = np.block(
        [[eigval * np.eye(A.shape[0]) - (D + t * (A - D)), x], [x.conj().T, 0]]
    )
    rhs = np.block([[(A - D) @ x], [0]])

    # Solve system_matrix * dy/dt = rhs for dy/dt
    dy_dt = np.linalg.solve(system_matrix, rhs)

    return dy_dt.T[0]


def main():
    n = 5

    # Define the matrices A and D
    A = np.random.rand(n, n)
    D = np.zeros((n, n))
    np.fill_diagonal(D, np.random.randn(n))

    # Define the initial conditions for each i
    y0_list = []
    for i in range(n):
        # Get the ith standard unit vector and diagonal entry of D
        e_i = np.zeros(n)
        e_i[i] = 1
        d_i = D[i, i]
        # Concatenate x(0) and eigval(0) as a real vector
        y0 = np.block([e_i, d_i])
        y0_list.append(y0)

    # Define the time span to solve the ODEs
    t_span = (0, 1)

    # Solve the ODEs for each i using solve_ivp
    solutions = []
    for y0 in y0_list:
        solution = spi.solve_ivp(f, t_span, y0, t_eval=[1], args=(A, D))
        solutions.append(solution)

    # Print the solutions for each i
    for i, solution in enumerate(solutions):
        print(f"Solution for i={i+1}:")
        print('eigval: ', solution.y[-1, -1])

    print(np.linalg.eig(A)[0])


if __name__ == "__main__":
    main()
