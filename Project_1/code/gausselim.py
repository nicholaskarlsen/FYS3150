# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
from __future__ import division  # Making sure integer division doenst sneak up on me
from numba import jit
import numpy as np
import time


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


@jit(nopython=True)
def general(n, a, b, c, timing=False):
    """"
    Solves for length n vector u in system Au = f(x)
    where A, [nxn] matrix and f(x) lengtth n vector
    and a function of x.

    takes input n, a, b, c
    n : number of steps
    a, b, c : Diagonal entries of sparse matrix A,
              below diagonal, diagonal and above
              diagonal respectively.
    """
    x = np.linspace(0, 1.0, n)
    h = 1.0 / (n - 1)      # Step size

    f = 100.0 * np.exp(-10.0 * x)  # Using directly to optimize for numba (jit)
    u = np.zeros(n)  # Initializing with np.zeros also sets Dirichlet bounds

    # Forward substitution
    for i in xrange(2, n - 1):  # 1 -> n-1
        b[i] -= ((a[i] * c[i - 1]) / b[i - 1])  # 3 * (n - 3)
        f[i] -= ((a[i] * f[i - 1]) / b[i - 1])  # 3 * (n - 3)
    u[-2] = f[-2] / b[-2]                       # 1
    # Backward substitution
    for i in xrange(n - 2, 0, -1):  # n-2 ->1
        u[i] = (f[i] - c[i] * u[i + 1]) / b[i]  # 3 * (n - 2)

    u *= h**2
    return x, u


@jit(nopython=True)
def specialized(n):
    """"
    Solves for length n vector u in system Au = f(x)
    where A, [nxn] matrix and f(x) lengtth n vector
    and a function of x.

    takes input n, a, b, c
    n : number of steps
    a, b, c : Diagonal entries of sparse matrix A,
              below diagonal, diagonal and above
              diagonal respectively.
    """

    x = np.linspace(0, 1.0, n)
    h = 1.0 / (n - 1)      # Step size

    f = 100.0 * np.exp(-10.0 * x)  # Using directly to optimize for numba (jit)
    u = np.zeros(n)  # Initializing with np.zeros also sets Dirichlet bounds
    b = np.ones(n) * 2

    # Forward substitution
    for i in xrange(2, n - 1):  # 1 -> n-1
        b[i] -= 1 / b[i - 1]            # 2 * (n-3)
        f[i] += f[i - 1] / b[i - 1]     # 2 * (n-3)
    u[-2] = f[-2] / b[-2]               # 1
    # Backward substitution
    for i in xrange(n - 2, 0, -1):  # n-2 ->1
        u[i] = (f[i] + u[i + 1]) / b[i]  # 2 * (n-2)

    u *= h**2

    return x, u


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    num = int(11)
    #x, u = specialized(num)

    a = np.ones(num) * -1      # Below diagonal
    b = np.ones(num) * 2       # Diagonal entries
    c = np.ones(num) * -1      # Above diagonal

    x, u = specialized(num)
    x2, u2 = general(num, a, b, c)
    print "general"
    print "x         numeric   exact"
    for i in xrange(num):
        print "%.4f" % x2[i], "--", "%.8f" % u2[i]

    print "Specialized:"
    print "x         numeric   exact"
    for i in xrange(num):
        print "%.4f" % x[i], "--", "%.8f" % u[i]
