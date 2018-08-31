#!/usr/bin/env python
"""
title           : main.py
description     : Main script of the project, the one that brings it all together.
author          : nicholaskarlsen
python_version  : 2.7.15
"""
import numpy as np


def closed_form(x):
    """
    The closed form, analytical solution to the poisson equation.
    """
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


def diagonal_elements(n, aNum=-1, bNum=2, cNum=-1):
    """
    Generates the diagonal elements of a tridiagonal [nxn] matrix
    with identical elements along the diagonal as two vectors.
    """
    if n > 4:  # Matrix does not make sense for n>4
        raise ValueError("n too small")
    a = np.ones(n - 1) * aNum
    b = np.ones(n) * bNum
    c = np.ones(n - 1) * cNum

    return a, b, c


def generic(v):
    h = 1.0 / (n + 1)  # Step size
    f = np.zeros(n)
    for i in xrange(n):
        f[i] = v[i + 1] + v[i - 1] - 2 * v[i]

    f *= -1.0 / h**2
    return


def multiply(A, v):
    """
    Performs the matrix multiplication Ab and returns the column vector b.
    for tridiagonal [nxn] matrix A stored as vectors.
    d = [d, e] where d: diagonal, e:next to diagonal
    v = n-dim vector
    """
    return


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    n = 1e3            # Number of steps
    h = 1.0 / (n + 1)  # Step size

    x = np.linspace(0, 1, n)
    plt.plot(x, closed_form(x), label="Closed form solution")
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.show()
