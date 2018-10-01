# Python Version: Python 2.7.15

from __future__ import division  # Nobody expects the integer division

import numpy as np
import matplotlib.pyplot as plt
import sys
from numba import jit
import main


def jacobi_rot(A, tol=1e-8, keep_counter=False):
    """
    Finds set of eigenvalues and its corresponding eigenvecotrs
    """

    print "--- Starting Jacobi method ---"
    counter = 0

    v = np.identity(len(A))    # Identidy matrix, used as set of orthonormal eigenvectors

    k, l = main.max_nondiag(A)

    while A[k, l]**2 >= tol or counter < len(A)**2:
        A, v = rotate(A, v, k, l)
        k, l = main.max_nondiag(A)
        counter += 1

    print "--- Found eigenvalues after %i rotations ---" % counter

    # in case i wish to store the counter for some reason
    if keep_counter is True:
        return A, v, counter
    else:
        return A, v


@jit(nopython=True)
def rotate(A, v, k, l):
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

    A[k, k] = c * c * a_kk - 2.0 * c * s * A[k, l] + s * s * a_ll
    A[l, l] = s * s * a_kk + 2.0 * c * s * A[k, l] + c * c * a_ll
    A[k, l] = 0.0  # hard-coding of the zeros
    A[l, k] = 0.0
    # and then we change the remaining elements
    for i in xrange(len(A)):
        if i != k and i != l:
            a_ik = A[i, k]
            a_il = A[i, l]
            A[i, k] = c * a_ik - s * a_il
            A[k, i] = A[i, k]
            A[i, l] = c * a_il + s * a_ik
            A[l, i] = A[i, l]
        # Compute new eigenvectors
        v_ik = v[i, k]
        v_il = v[i, l]
        v[i, k] = c * v_ik - s * v_il
        v[i, l] = c * v_il + s * v_ik

    return A, v


def testFunc():
    "Performs a series of tests to ensure proper functionality"
    N = 5
    eigTol = 1e-8

    X = main.construct(N)
    Y, v = jacobi_rot(X)
    jacobi_eigenvals = np.zeros(N)

    # Put the eigenvalues from the diagonal into a 1-D array
    for i in xrange(N):
        jacobi_eigenvals[i] = Y[i, i]

    # Sort list of eigenvalues & eigenvecs in ascending order
    permute1 = jacobi_eigenvals.argsort()
    jacobi_eigenvals = jacobi_eigenvals[permute1]
    jacobi_eigenvecs = v[:, permute1]

    Eigval_numpy, Eigvec_numpy = np.linalg.eig(X)
    permute2 = Eigval_numpy.argsort()
    Eigval_numpy = Eigval_numpy[permute2]
    Eigvec_numpy = Eigvec_numpy[:, permute2]
    # Checking that the computed eigenvalues match up with analytical eigenvalues
    for i in xrange(N):
        if abs(jacobi_eigenvals[i] - Eigval_numpy[i]) >= eigTol:
            print "-- Computed eigenvalues differ too much from numpy result, exiting --"
            sys.exit()
        else:
            pass

    print "Test function passed with tolerance %.1e." % eigTol
    print np.matmul(np.transpose(jacobi_eigenvecs[:,1]), jacobi_eigenvecs[:, 0])
    print np.matmul(np.transpose(Eigvec_numpy[:, 1]), Eigvec_numpy[:, 0])
    print 
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
