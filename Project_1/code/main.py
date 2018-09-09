from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import timeit
from gausselim import *
from lu import *
import time


def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


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
    if show == False:
        plt.close()
    else:
        plt.show()
    return


def analyticSolution(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


def relError(u, v):
    "Computes relative error between two sets of data u, v and returns max"
    err = np.abs((u[1:-2] - v[1:-2]) / v[1:-2])
    maxErr = np.amax(err)
    return maxErr


def ex_c(showplots=False):
    "Contains the calls pertaining to exercise 1c"

    List_of_N = [1e1, 1e2, 1e3, 1e4]
    plt.figure(figsize=(3.8, 3.8))
    for N in List_of_N:
        N = int(N)  # Float input causes errors
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1       # Above diagonal

        x_gauss, u_gauss = gauss_general(N, a, b, c)

        plt.plot(x_gauss, u_gauss, label="N = %i" % N)

    x_ana = np.linspace(0, 1, 1e3)  # lots of points to create smooth line
    u_ana = analyticSolution(x_ana)
    plt.plot(x_ana, u_ana, label="Analytic Solution")

    # plt.legend()
    # plt.savefig("../figs/ec1c_compare.png")

    figsetup(title="Testing for convergence of algorithm", xlab="x", ylab="u(x)", fname="ex1c_compare", show=True)
    return


def ex_d( showplots=False):

    general_times = []
    special_times = []

    List_of_N = [10, 25, 50, 75, 100, 250, 500, 750, 1e3, 2.5e3, 5e3, 7.5e3,
                 1e4, 1e5, 5e5, 1e6, 1e7, 1e8]

    for N in List_of_N:
        N = int(N)               # Float input causes errors
        # Preparing input vectors before function call
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1      # Above diagonal

        print "\nStarting timing for N=%.1E" % N

        """
        # wrapping function calls for timeit.timeit to work.
        def generalcall():
            gauss_general(N, a, b, c)

        def specialcall():
            gauss_specialized(N)

        print "Timing general algorithm"
        gen = timeit.timeit(generalcall, number=ncalls)
        general_times.append(gen)
        print "N=%.1E -> %.3E s per call" % (N, gen)
        print "Timing specialized algorithm"
        spe = timeit.timeit(specialcall, number=ncalls)
        special_times.append(spe)
        print "N=%.1E -> %.3E s per call" % (N, spe)
        """

        if N < 1E4:
            ncalls = 10000
        else:
            ncalls = 100


        t = 0
        for i in xrange(ncalls):
            t0 = time.time()
            gauss_general(N, a, b, c)
            t1 = time.time()
            t += t1 - t0
        general_times.append(t / ncalls)

        t = 0
        for i in xrange(ncalls):
            t0 = time.time()
            gauss_specialized(N)
            t1 = time.time()
            t += t1 - t0
        special_times.append(t / ncalls)

    log10N = np.log10(np.array(List_of_N))
    plt.figure(figsize=(3.8, 3.8))
    plt.plot(log10N, general_times, "x--", label="General algorithm")
    plt.plot(log10N, special_times, "x--", label="Specialized algorithm")
    figsetup(title="Average execution time of functions", xlab="$log_{10}N$",
             ylab="Time [s]", fname="ex1d_time", show=True)
    # Percentage difference in timings
    general_times = np.array(general_times)  # Convert to arrays for easier manipulation
    special_times = np.array(special_times)

    percent_diff = np.abs(general_times - special_times) / ((special_times + general_times) / 2.0) * 100

    plt.figure(figsize=(3.8, 3.8))
    plt.plot(log10N, percent_diff, "x--")
    figsetup(title="Percentage difference of execution time", xlab="$log_{10}N$",
             ylab="Percentage difference [%]", fname="ex1d_timediff", show=True)

    return


if __name__ == '__main__':
    #ex_c()
    ex_d()
