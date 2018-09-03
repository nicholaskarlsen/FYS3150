# An implementation of the Thomas Algorithm for solving
# Systems of the shape Ax = d, A=[n x n] tridiagonal matrix
# and x, d vectors.
from __future__ import division
import numpy as np


def thomas(a, b, c, n):

    for vec1 in [a, b, c]:
        for vec2 in [a, b, c]:
            if len(vec1) != len(vec2):
                raise ValueError("Input vectors of different length")

    c_dash = np.zeros(n)

    c_dash[0] = c[0] / b[0]

    for i in xrange(1, n):
        c_dash[i] = c[i] / (b[i] - a[i] * c_dash[i - 1])

    d_dash = np.zeros(n)

    d_dash[0] = d[0] / b[0]

    for i in xrange(1, n):
        d_dash[i] = (d[i] - a[i] * d_dash[i - 1]) / (b[i] - a[i] * c_dash[i - 1])

    x[-1] = d_dash[-1]

    for i in xrange(n - 1, 0, -1):  # range function does not include last entry
        x[i] = d_dash[i] - c_dash * x[i + 1]

    return

