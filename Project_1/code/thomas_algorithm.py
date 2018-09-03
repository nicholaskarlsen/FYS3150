# An implementation of the Thomas Algorithm for solving
# Systems of the shape Ax = d, A=[n x n] tridiagonal matrix
# and x, d vectors.
from __future__ import division
import numpy as np


def f(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


def thomas(a, b, c, x0=0, xn=1):

    n = len(b)

    x = np.zeros(n)
    d = np.zeros(n)
    c_dash = np.zeros(n)
    d_dash = np.zeros(n)

    x[0] = x0
    x[-1] = xn

    d[0] = f(x[0])

    c_dash[0] = c[0] / b[0]

    for i in xrange(1, n):  # Forward
        c_dash[i] = c[i] / (b[i] - a[i] * c_dash[i - 1])

    d_dash[0] = d[0] / b[0]

    for i in xrange(1, n):  # Backward
        d_dash[i] = (d[i] - a[i] * d_dash[i - 1]) / \
            (b[i] - a[i] * c_dash[i - 1])

    x[-1] = d_dash[-1]

    for i in xrange(n-2, 0, -1):
        x[i] = d_dash[i] - c_dash[i] * x[i + 1]

    return x


n = 10
a = np.ones(n) * -1
c = np.copy(a)
b = np.ones(n) * 2
