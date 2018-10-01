# Python Version: Python 2.7.15

from __future__ import division  # Nobody expects the integer division

import numpy as np
from numba import jit


@jit(nopython=True)
def max_nondiag(input_matrix):
    """
    Returns index i, j for the nondiagonal element which square is maximum.
    If multiple, equal maxima, return the last one
    Assumes symetrical matrix
    """
    N_input = len(input_matrix)
    currMax = 0  # Stores the current maximum
    for i in xrange(N_input):
        for j in xrange(i + 1, N_input):
            if input_matrix[i, j]**2 >= currMax:
                currMax = input_matrix[i, j]**2
                imax = i
                jmax = j
            else:
                pass

    return imax, jmax


@jit(nopython=True)
def rotate(mat, vec, k, l):
    "Performs Givens rotation"
    n_mat = len(mat)
    tau = (mat[l, l] - mat[k, k]) / (2.0 * mat[k, l])

    if tau >= 0:                                # tan(ang)
        t = 1.0 / (tau + np.sqrt(1.0 + tau**2))
    else:
        t = -1.0 / (-tau + np.sqrt(1.0 + tau**2))

    c = 1.0 / np.sqrt(1 + t**2)                 # cos(ang)
    s = t * c                                   # sin(ang)
    # update elements with indices k, l
    mat_kk = mat[k][k]
    mat_ll = mat[l][l]

    mat[k][k] = c * c * mat_kk - 2.0 * c * s * mat[k][l] + s * s * mat_ll
    mat[l][l] = s * s * mat_kk + 2.0 * c * s * mat[k][l] + c * c * mat_ll
    mat[k][l] = 0.0
    mat[l][k] = 0.0
    # update remaining elements
    for i in xrange(n_mat):
        # update eigenvector
        vec_ik = vec[i, k]
        vec_il = vec[i, l]
        vec[i, k] = c * vec_ik - s * vec_il
        vec[i, l] = c * vec_il + s * vec_ik
        if i != k and i != l:
            mat_ik = mat[i][k]
            mat_il = mat[i][l]
            mat[i][k] = c * mat_ik - s * mat_il
            mat[k][i] = mat[i][k]
            mat[i][l] = c * mat_il + s * mat_ik
            mat[l][i] = mat[i][l]

    return mat, vec


def jacobi_solve(matrix_input, tol=1e-8):
    # Create a copy of input matrix to avoid the original being changed outside of func.
    matrix = np.copy(matrix_input)
    dim = len(matrix)
    counter = 0     # Keeps track of no. transforms
    eigVec = np.identity(dim)  # used to form the eigenvectors

    p, q = max_nondiag(matrix)  # Fetch starting index

    while matrix[p, q]**2 >= tol:
        matrix, eigVec = rotate(mat=matrix, vec=eigVec, k=p, l=q)  # perform transform
        p, q = max_nondiag(matrix)  # Find indices of new max
        counter += 1
    print "Found eigenpair after %i transformations" % counter

    eigVal = np.zeros(dim)
    for i in xrange(dim):
        eigVal[i] = matrix[i, i]

    return eigVal, eigVec
