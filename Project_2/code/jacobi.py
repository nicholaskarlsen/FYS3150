# Python Version: Python 2.7.15

from __future__ import division  # Nobody expects the integer division

import numpy as np
import matplotlib.pyplot as plt
import sys
from numba import jit
import main


def jacobi_rot(A, tol=1e-8, keep_counter=False):
    """
    tol = Tolerance
    """

    print "--- Starting Jacobi method ---"
    counter = 0

    k, l = main.max_nondiag(A)

    while A[k, l]**2 >= tol or counter < len(A)**2:
        A = rotate(A, k, l)
        k, l = main.max_nondiag(A)
        counter += 1
        print counter

    print "--- Found eigenvalues after %i rotations ---" % counter

    if keep_counter is True:
        return A, counter
    else:
        return A


@jit(nopython=True)
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
    N = 5
    eigTol = 1e-8

    X = main.construct(N)
    Y = jacobi_rot(X)
    jacobi_eigenvals = np.zeros(N)

    # Put the eigenvalues from the diagonal into a 1-D array
    for i in xrange(N):
        jacobi_eigenvals[i] = Y[i, i]

    # Sort list of eigenvalues in ascending order
    jacobi_eigenvals = np.sort(jacobi_eigenvals)
    ana_eigenvals = np.sort(main.analyticalSolution(N))

    # Checking that the computed eigenvalues match up with analytical eigenvalues
    for i in xrange(N):
        if abs(jacobi_eigenvals[i] - ana_eigenvals[i]) >= eigTol:
            print "-- Computed eigenvalues differ too much from numpy result, exiting --"
            sys.exit()
        else:
            pass

    print "Test function passed with tolerance %.1e." % eigTol
    return


def num_rotations():
    N = np.linspace(5, 200, 10, dtype=int)
    no_rots = np.zeros(len(N))
    for i in range(len(N)):
        temp, no_rots[i] = jacobi_rot(main.construct(N[i]), keep_counter=True)

    a = np.polyfit(N, no_rots, deg=2)
    print a
    x = np.linspace(5, 1000, 1e5)
    extrapolation = a[0] * x**2 + a[1] * x + a[2]

    plt.figure(figsize=(3.8, 3.8))
    plt.loglog(N, no_rots, "x", label="Data")
    plt.loglog(x, extrapolation, label="Etrapolation")

    main.figsetup(title="Number of itterations needed for N-dim",
                  xlab="N", ylab="No. Rotations", fname="norots")

    return


if __name__ == '__main__':
    testFunc()
    # num_rotations()

    """

    A_test = main.construct(50)

    B = jacobi_rot(A_test)

    eigenvals = np.zeros(len(B))
    for i in range(len(B)):
        eigenvals[i] = B[i, i]

    print np.sort(np.linalg.eigvals(B))
    print np.sort(eigenvals)
    
    """
