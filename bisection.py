# Reference: Numerical Methods for Eigenvalue Problems (2012) Ch. 6
# Additional reference for sturm sequences: https://www.cse.psu.edu/~b58/cse456/lecture13.pdf

# matrix must be real symmetric tridiagonal

import numpy as np


def sturm_evaluate(t, a, b):
    p = [1, t - a[0]]
    c = 0
    if p[0] * p[1] < 0 or p[1] == 0:
        c += 1
    n = len(a)
    for m in range(2, n+1):
        p.append((t - a[m-1]) * p[m-1] - abs(b[m-2]) ** 2 * p[m-2])
        if p[m] * p[m-1] < 0 or p[m] == 0:
            c += 1
    return c


def sturm_bisection(k, a, b, alpha, beta):
    n = a.shape[0]
    epsilon = 1e-6
    while beta - alpha > epsilon:
        gamma = (beta + alpha) / 2
        c = sturm_evaluate(gamma, a, b)
        print(k, c)

        if k <= n - c:
            beta = gamma
        else:
            alpha = gamma
    print()
    return gamma


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
    T = symmetric_tridiagonal_matrix(10)
    a = np.diag(T)
    b = np.diag(T, -1)

    alpha, beta = gershgorin_bound(a, b)

    num_eigs = T.shape[0]
    eigs = []

    for k in range(1, num_eigs+1):
        eig = sturm_bisection(k, a, b, alpha, beta)
        eigs.append(eig)
        print()

    print(eigs)
    print(sorted(np.linalg.eig(T)[0]))


main()