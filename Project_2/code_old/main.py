# Generates all of the figures and data for the report
# Python Version: Python 2.7.15

from __future__ import division  # Nobody expects the integer division

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.sparse
from numba import jit
import jacobi


def figsetup(title, xlab, ylab, fname, show=False):
    """
    Sets up and saves figure for usage in report
    usage:
    plot(...)
    plot(...)
    figsetup("filename")
    """
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(fname)
    plt.tight_layout()
    plt.title(title)
    plt.legend()
    plt.savefig("../figs/" + fname + ".png", dpi=250)
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def construct(N, rhoMax=0, potential=False):
    """
    Generates the matrix A as defined in report.
    """
    if potential == 1:
        i = np.linspace(0, rhoMax, N)
        h = (rhoMax) / (N + 1)
        rho = i * h
        d = 2.0 / h**2 * np.ones(N) + rho**2

    if potential is False:
        h = 1.0 / (N + 1)
        d = 2.0 / h**2 * np.ones(N) + potential
    a = - 1.0 / h**2 * np.ones(N)
    A = scipy.sparse.spdiags([a, d, a], [-1, 0, 1], N, N).toarray()  # Generates matrix
    return A


@jit
def max_nondiag(input_matrix, tol=1e-8):
    """
    Returns index i, j for the nondiagonal element which square is maximum.
    If multiple, equal maxima, return the last one
    Assumes symetrical matrix
    """
    temp_matrix = np.copy(input_matrix)**2

    currMax = 0  # Stores the current maximum
    for i in xrange(len(input_matrix)):
        for j in xrange(i + 1, len(input_matrix)):
            if temp_matrix[i, j]**2 >= currMax:
                currMax = temp_matrix[i, j]**2
                imax = i
                jmax = j
            else:
                pass

    return imax, jmax


def analyticalSolution(N):

    return lambd


def exercise_2b():
    """ Contains calls relevant to question 2b """

    return


def test_max_dondiag():
    testMat = np.array([[10, 2, 4, 6],
                        [2, 11, -9, 5],
                        [4, -9, 12, 3],
                        [6, 5, 3, 13], ])

    i, j = max_nondiag(testMat)

    if testMat[i, j]**2 != 81:
        raise "Found wrong entry"

    if i == j:
        raise "Found diagonal entry"

    return


def ex_d():
    "contains the calls pertaining to exercise 2e"
    n = 400
    X = construct(n, rhoMax=1, potential=1)
    Y, Eigvec_jacobi = jacobi.jacobi_rot(X)

    Eigval_jacobi = np.zeros(n)
    for i in xrange(n):
        Eigval_jacobi[i] = Y[i, i]
    Eigval_numpy, Eigvec_numpy = np.linalg.eig(X)
    permute = Eigval_numpy.argsort()
    Eigval_numpy = Eigval_numpy[permute]
    Eigvec_numpy = Eigvec_numpy[:, permute]

    permute = Eigval_jacobi.argsort()
    Eigval_jacobi = Eigval_jacobi[permute]
    Eigvec_jacobi = Eigvec_jacobi[:, permute]
    plt.plot(Eigvec_numpy[1])
    plt.show()
    return


if __name__ == '__main__':
    # test_max_dondiag()
    ex_d()
