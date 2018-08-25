import numpy as np
# import matplotlib.pyplot as plt

# Initial Conditions & parameters
n = 1e3
h = 1.0 / float(n)
u0 = 0  # u(0) = 0
u1 = 0  # u(1) = 0
x = np.linspace(0, 1, h)


def tridiagonal(n, dNum=2, eNum=-1):
    """
    Generates the diagonal elements of a tridiagonal [nxn] matrix
    with identical elements along the diagonal as two vectors.
    """
    if n < 4:  # Matrix does not make sense for n>4
        raise ValueError("n too small")

    d = np.ones(n) * dNum
    e = np.ones(n - 2) * eNum

    return [d, e]


def multiply(A, v):
    """
    Performs the matrix multiplication Ab and returns the column vector b.
    for tridiagonal [nxn] matrix A stored as vectors.
    d = [d, e] where d: diagonal, e:next to diagonal
    v = n-dim vector
    """

    if np.shape(A)[0] != len(v):
        raise ValueError("Can not multiply due to dimensions.")

    n = len(v)
    m = np.shape(A)[1]

    b = np.zeros(n)

    for i in xrange(n):
        for j in xrange(m):
            b[i] += A[i][j] * v[i]

    return b


A = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1],
              ])

v = np.array([1, 2, 3, 4])

print multiply(A, v)
