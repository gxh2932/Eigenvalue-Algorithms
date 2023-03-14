#Reference: Numerical Methods for Eigenvalue Problems (2012) Ch. 6

# matrix must be real symmetric tridiagonal

import numpy as np


def construct_characteristic_polynomial(a, b):
    '''
    Constructs the characteristic polynomial of a real symmetric tridiagonal matrix. Uses the recurrence relation
    p_n(x) = (a_n - x)p_{n-1}(x) - b_{n-1}^2 p_{n-2}(x).
    :param a:
    :param b:
    :return:
    '''
    n = a.shape[0]
    p = [1, np.poly1d([-1, a[0]])]
    for i in range(1, n):
        p_1 = np.polymul(p[i], [-1, a[i]])
        p_2 = np.polymul(p[i-1], [-b[i-1]**2])
        p.append(np.polyadd(p_1, p_2))
    return p[-1]


def laguerre_method(p, x0, epsilon=1e-6, max_iter=100):
    """
    Implements the Laguerre method for finding a root of a polynomial function.

    Parameters:
    p (np.ndarray): Coefficients of the polynomial function, in decreasing order.
    dp (np.ndarray): Coefficients of the first derivative of the polynomial function, in decreasing order.
    d2p (np.ndarray): Coefficients of the second derivative of the polynomial function, in decreasing order.
    x0 (float): Initial guess for the root.
    epsilon (float): Tolerance for convergence.
    max_iter (int): Maximum number of iterations.

    Returns:
    float: Approximation of the root.
    """
    x = x0
    n = len(p)

    dp = np.polyder(p)
    d2p = np.polyder(dp)

    for i in range(max_iter):
        p_val = np.polyval(p, x)
        dp_val = np.polyval(dp, x)
        d2p_val = np.polyval(d2p, x)

        if abs(p_val) < epsilon:
            return x

        G = dp_val / p_val
        H = G**2 - d2p_val / p_val

        if G >= 0:
            a = n / (G + np.emath.sqrt((n - 1) * (n * H - G**2)))
        else:
            a = n / (G - np.emath.sqrt((n - 1) * (n * H - G**2)))

        x -= a

        if abs(a) < epsilon:
            return x

    Exception('Maximum number of iterations exceeded.')


def gershgorin_bound(a, b):
    n = a.shape[0]
    alpha = np.min([a[0]-np.abs(b[0]), a[n-1]-np.abs(b[n-2])] + [a[i]-np.abs(b[i])-np.abs(b[i-1]) for i in range(1, n-1)])
    beta = np.max([a[0]+np.abs(b[0]), a[n-1]+np.abs(b[n-2])] + [a[i]+np.abs(b[i])+np.abs(b[i-1]) for i in range(1, n-1)])
    return alpha, beta


def symmetric_tridiagonal_matrix(n):
    d = np.random.rand(n)
    e = np.random.rand(n-1)
    T = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
    return T


def main():
    T = symmetric_tridiagonal_matrix(30)
    a = np.diag(T)
    b = np.diag(T, -1)

    alpha, beta = gershgorin_bound(a, b)

    num_eigs = T.shape[0]
    eigs = []

    p = construct_characteristic_polynomial(a, b)
    x_0 = (alpha + beta) / 2

    for k in range(1, num_eigs+1):
        eig = laguerre_method(p, x_0)
        eigs.append(eig)

        # update the characteristic polynomial
        p = np.polydiv(p, [-1, eig])[0]

    print(sorted(eigs))
    print(sorted(np.linalg.eig(T)[0]))


main()