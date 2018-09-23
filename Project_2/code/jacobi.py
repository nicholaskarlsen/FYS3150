# Python Version: Python 2.7.15
from __future__ import division  # Defaults to float division
import numpy as np
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


if __name__ == '__main__':

    d = np.ones(N - 2) * 2.0 / h**2
    a = - np.ones(N - 2) * 1.0 / h**2

    A_sparse = scipy.sparse.spdiags([a, d, a], [-1, 0, 1], N_test - 2, N_test - 2).toarray()  # Generates matrix

    print jacobi_rot(A=A_test, N=N_test)
