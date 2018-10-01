from __future__ import division  # Nobody expects the integer division

import numpy as np
import jacobi_eigensolver
import main


def test_max_dondiag():
    testMat = np.array([[10, 2, 4, 6],
                        [2, 11, -9, 5],
                        [4, -9, 12, 3],
                        [6, 5, 3, 13], ])

    i, j = jacobi_eigensolver.max_nondiag(testMat)

    if testMat[i, j]**2 != 81:
        raise ValueError("Found incorrect entry")

    if i == j:
        raise ValueError("Found diagonal entry")

    return


def test_construct_matrix():
    "If i construct the correct matrix, numpy should yield analytical eigenvals"
    N = 10
    A, rho = main.construct_matrix(N)

    eval_np, evec_np = np.linalg.eig(A)  # Calculate eigenpair using numpy
    eval_np, evec_np = main.sort_eigenpair(eval_np, evec_np)  # sorts to ascending order

    # Get analytical eigenvalues
    h = 1.0 / (N + 1)
    d = 2.0 / h**2
    a = 1.0 / h**2
    j = np.linspace(1, N, N)  # j = 1, 2, ... , N-1
    lambd = np.sort(d + 2 * a * np.cos(j * np.pi / (N + 1)))

    for i in xrange(N):
        if abs(eval_np[i] - lambd[i]) > 1e-8:
            raise ValueError("Eigenvalues of constructed matrix differ from analytical")

    return


def test_jacobi_solve():
    """
    If jacobi_solve() works, i expect:
    (1) Get the same eigenvalues as numpy (& by exctension, analytical. see construct_matrix test)
    (2) The resultant eigenvectors must be orhogonal.
    """
    N = 6
    A, rho = main.construct_matrix(N)  # Generate matrix

    eval_np, evec_np = np.linalg.eig(A)  # Calculate eigenpair using numpy
    eval_np, evec_np = main.sort_eigenpair(eval_np, evec_np)  # sorts to ascending order

    eval_ja, evec_ja = jacobi_eigensolver.jacobi_solve(A)  # Calculate egeinpair using jacobi method
    eval_ja, evec_ja = main.sort_eigenpair(eval_ja, evec_ja)  # sort to ascending order

    # Check that computed eigenvalues correspond to numpy equivalent
    for i in xrange(N):
        if abs(eval_ja[i] - eval_np[i]) > 1e-8:
            raise ValueError("Computed eigenvalues differ too much from reference")

    # Check orthogonality of eigenvectors
    for i in xrange(N):
        for j in xrange(N):
            inner_prod = np.matmul(np.transpose(evec_ja[:, i]), evec_ja[:, j])
            # checking that inner product for i, j => kronicker delta within tolerance
            if i != j and inner_prod > 1e-8:
                raise ValueError("orthogonality of eigenvectors not retained")
            if i == j and abs(inner_prod - 1) > 1e-8:
                raise ValueError("orthogonality of eigenvectors not retained")
    return


def run_tests():
    "Runs all tests in sequence"
    test_max_dondiag()
    test_construct_matrix()
    test_jacobi_solve()

    print "All tests completed."
    return


if __name__ == '__main__':
    run_tests()
