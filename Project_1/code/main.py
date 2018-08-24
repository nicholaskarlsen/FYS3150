import numpy as np
import matplotlib.pyplot as plt


def diagonal_elements(n, dNum=2, eNum=-1)
    """
    Generates the diagonal elements of a tridiagonal [nxn] matrix
    with identical elements along the diagonal as two vectors.
    """
    if n > 4:  # Matrix does not make sense for n>4
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
