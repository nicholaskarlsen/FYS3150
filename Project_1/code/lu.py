# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
# Performs LU decomposition to solve tridiagonal matrix using scipy
import scipy.sparse
import scipy.linalg
import numpy as np
import main


def LU_benchmark(n):
    "Solves system by LU decomposition for a given n using scipy"
    n = int(n)  # In case input is float, needs to be int or else errors.
    h = 1.0 / float(n)  # Step size

    a = np.ones(n) * (-1)  # Entries below diagonal
    b = np.ones(n) * 2     # Diagonal entries
    c = np.copy(a)         # Above diagonal, copy of a.

    A = scipy.sparse.spdiags([a, b, c], [-1, 0, 1], n, n).toarray()  # Generates matrix

    xvals = np.linspace(0, 1, n)
    f = main.f_func(xvals)  # RHS of equation

    LU = scipy.linalg.lu_factor(A)  # LU decomposition
    v = scipy.linalg.lu_solve(LU, f)  # Solving system

    v *= h**2   # Revert scaling for solutiuon

    return xvals, v


if __name__ == '__main__':
    # Testbench ensuring everything works
    import time
    import matplotlib.pyplot as plt
    times = []
    nvals = [1e1, 1e2, 1e3, 1e4]
    for n in nvals:
        t0 = time.time()
        LU_benchmark(n)
        t1 = time.time()
        times.append(t1 - t0)

    plt.plot(np.log10(nvals), times, "x--", label="Scipy LU-Decomposition")
    plt.ylabel("Time [s]")
    plt.xlabel("$log_{10}n$")
    plt.legend()
    plt.show()
