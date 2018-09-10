# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
from __future__ import division  # Making sure integer division doenst sneak up on me
import numpy as np
from numba import jit
import sys
from main import *


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


@jit(nopython=True)
def gauss_general(n, a, b, c):
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

    x = np.linspace(0, 1.0, n)    # x in [0, 1]
    h = 1.0 / n      # Step size

    f = 100.0 * np.exp(-10.0 * x)  # Using directly to optimize for numba (jit)
    u = np.zeros(n)  # Initializing with np.zeros also sets Dirichlet bounds

    # Forward substitution
    for i in xrange(2, n-1):  # 1 -> n-1
        b[i] -= ((a[i] * c[i - 1]) / b[i - 1])
        f[i] -= ((a[i] * f[i - 1]) / b[i - 1])
    u[-2] = f[-2] / b[-2]
    # Backward substitution
    for i in xrange(n - 1, 0, -1):  # n-2 ->1
        u[i] = (f[i] - c[i] * u[i + 1]) / b[i]

    u *= h**2

    return x, u


@jit(nopython=True)
def gauss_specialized(n):

    x = np.linspace(0, 1, n)    # x in [0, 1]
    h = (x[-1] - x[0]) / n

    f = 100 * np.exp(-10 * x)  # Using directly to optimize
    u = np.zeros(n)


    "Tried to calculate _b outside of loop, resulted in lower performance"
    #i_array = np.linspace(1, n-1, n)
    #_b = (i_array + 1.0) / i_array

    # Forward
    for i in xrange(2, n):
        #_b[i] = (i + 1) / i  # Performing this outside of loop -> slower
        #_f[i] = f[i] + ((i - 1.0) * _f[i - 1]) / i
        b[i] = 2 - 1 / b[i - 1]
        f[i] = f[i] + f[i - 1] / b[i - 1]

    u[-2] = _f[-2] / _b[-2]
    # Backward
    for i in xrange(n - 1, 0, -1):
        u[i] = i / (i + 1) * (_f[i] + u[i + 1])

    u *= h**2

    return x, u


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    num = int(10)
    #x, u = gauss_specialized(num)

    a = np.ones(num) * -1      # Below diagonal
    b = np.ones(num) * 2       # Diagonal entries
    c = np.ones(num) * -1      # Above diagonal

    x2, u2 = gauss_general(num, a, b, c)
    ana = analyticSolution(x2)

    #plt.plot(x, u, label="Specialized")
    plt.plot(x2, u2, label="General")
    plt.plot(x2, ana, label="Analytic")
    plt.legend()
    plt.show()
