# Python Version: Python 2.7.15
from __future__ import division  # Defaults to float division
import numpy as np
import scipy.sparse
import main


def jacobi_rot(A, tol=1e-8):
    """
    tol = Tolerance
    """

    print "-- Starting Rotation --"

    counter = 0

    k, l = main.max_nondiag(A)

    while A[k, l]**2 >= tol:
        A = rotate(A, k, l)
        k, l = main.max_nondiag(A)
        counter += 1

    print "-- Rotation done --"
    print "Found eigenvalues after %i rotations" % counter
    return A


def rotate(A, k, l):
    """ Givens rotation """
    tau = (A[l, l] - A[k, k]) / (2 * A[k, l])   # cot(2*ang)

    if tau >= 0:                                # tan(ang)
        t = 1.0 / (tau + np.sqrt(1.0 + tau**2))
    else:
        t = -1.0 / (-tau + np.sqrt(1.0 + tau**2))

    c = 1.0 / np.sqrt(1 + t**2)                 # cos(ang)
    s = t * c                                   # sin(ang)

    a_kk = A[k, k]
    a_ll = A[l, l]

    A[k][k] = c * c * a_kk - 2.0 * c * s * A[k][l] + s * s * a_ll
    A[l][l] = s * s * a_kk + 2.0 * c * s * A[k][l] + c * c * a_ll
    A[k][l] = 0.0  # hard-coding of the zeros
    A[l][k] = 0.0
    # and then we change the remaining elements
    for i in xrange(len(A)):
        if i != k and i != l:
            a_ik = A[i][k]
            a_il = A[i][l]
            A[i][k] = c * a_ik - s * a_il
            A[k][i] = A[i][k]
            A[i][l] = c * a_il + s * a_ik
            A[l][i] = A[i][l]

    return A


def testFunc():
    "Performs a series of tests to ensure proper functionality"
    X = main.construct(10)
    Y = jacobi_rot(X)
    return


if __name__ == '__main__':
    testFunc()

    """

    A_test = main.construct(50)

    B = jacobi_rot(A_test)

    eigenvals = np.zeros(len(B))
    for i in range(len(B)):
        eigenvals[i] = B[i, i]

    print np.sort(np.linalg.eigvals(B))
    print np.sort(eigenvals)
    
    """
