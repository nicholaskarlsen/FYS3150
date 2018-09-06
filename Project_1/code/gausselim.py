# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
from __future__ import division
import numpy as np


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


def gauss_general(n, sourceterm, aVal=-1, bVal=2, cVal=-1):

    h = 1 / (n + 1)

    a = np.ones(n) * aVal       # Below diagonal
    b = np.ones(n) * bVal       # Diagonal entries
    c = np.ones(n) * cVal       # Above diagonal

    x = np.linspace(0, 1, n)    # x in [0, 1]

    f = sourceterm(x)

    _f = np.zeros(n)
    _b = np.zeros(n)
    u = np.zeros(n)

    _b[0] = b[0]
    _f[0] = f[0]

    # Forward step
    for i in range(2, n - 1):
        _b[i] = b[i] - (a[i] * c[i - 1]) / b[i - 1]
        _f[i] = f[i] - (a[i] * _f[i - 1]) / b[i - 1]
    # Backward Step
    for i in range(n - 1, 1, -1):
        u[i - 1] = (_f[i - 1] - c[i - 1] * u[i]) / _b[i - 1]

    u[-1] = _f[-1] / _b[-1]

    u *= h**2

    return x, u


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    x, u = gauss_general(10, f_func)

    plt.plot(x, u)
    plt.show()
