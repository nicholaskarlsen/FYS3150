# Generates all of the figures and data for the report
# Python Version: Python 2.7.15

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.sparse


def constructA(N):
    """
    Generates the matrix A as defined in report.
    """
    h = 1.0 / N
    d = 2.0 / h**2 * np.ones(N - 2)
    a = - 1.0 / h**2 * np.ones(N - 2)
    A = scipy.sparse.spdiags([a, d, a], [-1, 0, 1], N - 2, N - 2).toarray()  # Generates matrix

    eig = np.linalg.eigvals(A)

    return A, eig


def max_nondiag(A, tol=1e-8):
    " Finds the maximum of the squares of nondiagonal elements in A "
    A_copy = np.copy(A)             # Creates a copy of A
    np.fill_diagonal(A_copy, 0)     # Fills diagonal of A copy with 0

    return np.amax(A_copy**2)


def analyticalSolution(N):
    h = 1.0 / N
    d = 2.0 / h**2
    a = 1.0 / h**2
    j = np.linspace(1, N - 1, N - 1)  # j = 1, 2, ... , N-1
    lambd = d + 2 * a * np.cos(j * np.pi / (N + 1))
    return lambd


def librarySolution(N):
    """
    Finds eigenvalues using the diagonalization functions in numpy
    """
    A = constructA(N)
    A_diag = np.diag(A)
    return


def exercise_2b():
    """ Contains calls relevant to question 2b """

    return


if __name__ == '__main__':
    A, eig = constructA(10)
    eig2 = analyticalSolution(10)
    print eig[0]
    print eig2[0]
    print "Done"

    print max_nondiag(A)
