import scipy.sparse
import scipy.linalg
import numpy as np


def LU_benchmark(n):
    h = 1.0 / float(n)

    a = np.ones(n) * (-1)
    b = np.ones(n) * 2
    c = np.copy(a)

    A = scipy.sparse.spdiags([a, b, c], [-1, 0, 1], n, n).toarray()

    xvals = np.linspace(0, 1, n)

    def f_func(x):
        """ Source term """
        return 100 * np.exp(-10 * x)

    f = f_func(xvals)

    LU = scipy.linalg.lu_factor(A)

    v = scipy.linalg.lu_solve(LU, f)

    v *= h**2

    return xvals, v


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    times = []
    nvals = [10, 100, 10000, 10000]
    for n in nvals:
        t0 = time.time()
        LU_benchmark(n)
        t1 = time.time()
        times.append(t1 - t0)

    plt.plot(nvals, times)
    plt.show()
