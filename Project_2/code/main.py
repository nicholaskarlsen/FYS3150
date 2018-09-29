# Generates all of the figures and data for the report
# Python Version: Python 2.7.15

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.sparse


def construct(N):
    """
    Generates the matrix A as defined in report.
    """
    h = 1.0 / (N + 1)
    d = 2.0 / h**2 * np.ones(N)
    a = - 1.0 / h**2 * np.ones(N)
    A = scipy.sparse.spdiags([a, d, a], [-1, 0, 1], N, N).toarray()  # Generates matrix
    return A


def max_nondiag(input_matrix, tol=1e-8):
    """
    Returns index i, j for the nondiagonal element which square is maximum.
    If multiple, equal maxima, return the last one
    """
    temp_matrix = np.copy(input_matrix)  # Creates a copy to manipulate
    np.fill_diagonal(temp_matrix, 0)     # Fills diagonal of temp with 0

    currMax = 0  # Stores the current maximum
    for i in xrange(len(input_matrix)):
        for j in xrange(len(input_matrix)):
            if temp_matrix[i, j]**2 >= currMax:
                currMax = temp_matrix[i, j]**2
                imax = i
                jmax = j
            else:
                pass

    if imax == jmax:
        raise "max_nondiag() found a diagonal element"

    return imax, jmax


def analyticalSolution(N):
    h = 1.0 / N
    d = 2.0 / h**2
    a = 1.0 / h**2
    j = np.linspace(1, N - 1, N - 1)  # j = 1, 2, ... , N-1
    lambd = d + 2 * a * np.cos(j * np.pi / (N + 1))
    return lambd


def exercise_2b():
    """ Contains calls relevant to question 2b """

    return


def test_max_dondiag():
    testMat = np.array([[1, 9, 4, 5],
                        [1, 11, -9, 4],
                        [1, 3, 10, 8],
                        [1, 3, 4, 5]])

    i, j = max_nondiag(testMat)

    if i != 1 and j != 2:
        raise "max_nondiag() returned unexpected indices"
    return


if __name__ == '__main__':
    test_max_dondiag()
