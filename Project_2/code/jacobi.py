# Python Version: Python 2.7.15
from __future__ import division  # Defaults to float division
import numpy as np
import scipy.sparse
import main


def jacobi_rot(A, N, tol=1e-8):
    """
    tol = Tolerance
    """

    S = np.zeros([N, N], dtype=np.float64)

    def tau():
        return (a[l, l] - a[k, k]) / (2 * a[k, l])

    while main.max_nondiag(A) >= tol:
        pass

    return S


def rotate(A, k, l):
    sim_A = np.copy(A)  # Creates a copy of A to modify

    tau = (A[l, l] - A[k, k]) / (2 * A[k, l])   # cot(2*ang)

    t = - tau + np.sqrt(1 + tau**2)
    t = - tau - np.sqrt(1 + tau**2)

    c = 1.0 / np.sqrt(1 + t**2)                 # cos(ang)
    s = t * c                                   # sin(ang)

    sim_A[k, k] = A[k, k] * c**2 - 2 * A[k, l] * c * s + A[l, l] * s**2
    sim_A[l, l] = A[l, l] * c**2 + 2 * A[k, l] * c * s + A[k, k] * s**2
    sim_A[k, l] = (a[k, k] - a[l, l]) * c * s + a[k, l] * (c**2 - s**2)

    return


if __name__ == '__main__':

    N_test = 100

    h = 1.0 / N_test

    d = np.ones(N_test - 2) * 2.0 / h**2
    a = - np.ones(N_test - 2) * 1.0 / h**2

    A_test = scipy.sparse.spdiags(
        [a, d, a], [-1, 0, 1], N_test - 2, N_test - 2).toarray()  # Generates matrix
