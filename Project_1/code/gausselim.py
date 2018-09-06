# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
from __future__ import division
import numpy as np
from main import *


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


def gauss_general(n, sourceterm, aVal=-1, bVal=2, cVal=-1):

    a = np.ones(n, dtype=float) * aVal       # Below diagonal
    b = np.ones(n, dtype=float) * bVal       # Diagonal entries
    c = np.ones(n, dtype=float) * cVal       # Above diagonal

    x = np.linspace(0, 1, n)    # x in [0, 1]
    h = (x[-1] - x[0]) / n

    f = sourceterm(x)

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
        #u[i - 1] = (_f[i - 1] - c[i - 1] * u[i]) / _b[i]
        u[i] = (_f[i] - c[i] * u[i + 1]) / _b[i]

    u *= h**2

    return x, u


def gauss_specialized():

    return


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    num = int(1e6)
    x, u = gauss_general(num, f_func)
    ana = analyticSolution(x)

    plt.plot(x, u, label="Gauss elim")
    plt.plot(x, ana, label="Analytic")
    plt.legend()
    plt.show()
