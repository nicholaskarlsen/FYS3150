# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
from __future__ import division
import numpy as np
from numba import jit
from main import *
import timeit


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)

@jit(nopython=True)
def gauss_general(n, aVal=-1, bVal=2, cVal=-1):

    a = np.ones(n) * aVal       # Below diagonal
    b = np.ones(n) * bVal       # Diagonal entries
    c = np.ones(n) * cVal       # Above diagonal

    x = np.linspace(0, 1, n)    # x in [0, 1]
    h = (x[-1] - x[0]) / n

    # f = sourceterm(x)
    f = 100 * np.exp(-10 * x)  # Using directly to optimize

    _f = np.zeros(n)
    _b = np.zeros(n)
    u = np.zeros(n)  # Initializing also sets Dirichlet bounds

    # If different boundary conditions are desired
    # u[0] = ...
    # u[-1] = ...

    _b[0] = b[0]
    _f[0] = f[0]

    # Forward substitution
    for i in range(1, n):
        _b[i] = b[i] - (a[i - 1] * c[i - 1]) / _b[i - 1]
        _f[i] = f[i] - (a[i] * _f[i - 1]) / _b[i - 1]

    u[n - 1] = _f[n - 1] / _b[n - 1]

    # Backward substitution
    for i in range(n - 2, 0, -1):
        u[i] = (_f[i] - c[i] * u[i + 1]) / _b[i]

    u *= h**2

    return x, u

@jit(nopython=True)
def gauss_specialized(n):

    x = np.linspace(0, 1, n)    # x in [0, 1]
    h = (x[-1] - x[0]) / n

    # f = sourceterm(x)
    f = 100 * np.exp(-10 * x)  # Using directly to optimize

    _f = np.zeros(n)
    _b = np.zeros(n)
    u = np.zeros(n)

    _b[0] = b[0]
    _f[0] = f[0]

    # Forward
    for i in range(1, n):
        _b[i] = b[i] - (a[i - 1] * c[i - 1]) / _b[i - 1]
        _f[i] = f[i] - (a[i] * _f[i - 1]) / _b[i - 1]

    u[n - 1] = _f[n - 1] / _b[n - 1]
    # Backward
    for i in range(n - 2, 0, -1):
        u[i] = (_f[i] - c[i] * u[i + 1]) / _b[i]

    u *= h**2

    return x, u


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    num = int(1e6)
    x, u = gauss_general(num)
    ana = analyticSolution(x)

    plt.plot(x, u, label="Gauss elim")
    plt.plot(x, ana, label="Analytic")
    plt.legend()
    plt.close()
