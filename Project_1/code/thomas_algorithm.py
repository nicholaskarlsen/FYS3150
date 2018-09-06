# An implementation of the Thomas Algorithm for solving
# Systems of the shape Ax = d, A=[n x n] tridiagonal matrix
# and x, d vectors.
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


def thomas(n, x0=0, xn=1):

    h = 1 / (n + 1)

    a = np.ones(n) * -1
    c = np.copy(a)
    b = np.ones(n) * 2

    pos = np.linspace(0, 1, n)
    x = np.zeros(n)
    d = np.zeros(n)
    c_dash = np.zeros(n)
    d_dash = np.zeros(n)

    x[0] = x0
    x[-1] = xn

    d = f_func(pos)

    c_dash[0] = c[0] / b[0]

    for i in xrange(1, n):  # Forward
        c_dash[i] = c[i] / (b[i] - a[i] * c_dash[i - 1])

    d_dash[0] = d[0] / b[0]

    for i in xrange(1, n):  # Backward
        d_dash[i] = (d[i] - a[i] * d_dash[i - 1]) / \
            (b[i] - a[i] * c_dash[i - 1])

    x[-1] = d_dash[-1]

    for i in xrange(n - 2, 0, -1):
        x[i] = d_dash[i] - c_dash[i] * x[i + 1]

    return x * h**2, pos


def thomas_alt(n):

    xvals = np.linspace(0, 1, n)

    a = np.ones(n) * (-1)
    b = np.ones(n) * 2
    c = np.ones(n) * (-1)
    v = np.zeros(n)

    f = np.zeros(n)
    f_dash = np.zeros(n)

    f_dash[1] = f_func(0)

    for i in xrange(1, n):
        b[i] = b[i] - a[i] * c[i - 1] / b[i - 1]
        f_dash[i] = f[i] - a[i] * f_dash[i - 1] / b[i - 1]

    v[-1] = f_dash[-1] / b[-1]

    for i in xrange(n - 2, 0, -1):
        print(i)
        v[i] = (f_dash[i] - c[i] * v[i + 1]) / b[i]

    return v, xvals

if __name__ == '__main__':
    s = thomas_alt(100)

    plt.plot(s[1], s[0])
    plt.show()