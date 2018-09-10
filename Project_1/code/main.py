# By Nicholas Karlsen
# Python 2.7.14 (Anaconda)
"""
Contains all of the function calls to produce material for Project 1 as well as some
other useful functions.
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from gausselim import *
from lu import *


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
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def analyticSolution(x):
    return 1.0 - (1.0 - np.exp(-10.0)) * x - np.exp(-10.0 * x)


def relError(v, u):
    "Computes relative error between two sets of data u, v and returns max"
    # Dont include boundaries as they are equal and zero (-> error)
    err = abs((v[1:-1] - u[1:-1]) / u[1:-1])
    maxErr = np.amax(err)
    return maxErr


""" Below are functions containing calls to generate plots etc for report """


def ex_c(showplots=False):
    "Contains the calls pertaining to exercise 1c"

    print "starting (c)"
    List_of_N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    plt.figure(figsize=(3.8, 3.8))
    for N in List_of_N:
        print "N=%.1E" % N
        N = int(N)  # Float input causes errors
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1       # Above diagonal

        x_gauss, u_gauss = gauss_general(N, a, b, c)

        plt.plot(x_gauss, u_gauss, label="N = %.1E" % N)

    x_ana = np.linspace(0, 1, 1e3)  # lots of points to create smooth line
    u_ana = analyticSolution(x_ana)
    plt.plot(x_ana, u_ana, label="Analytic Solution")

    # plt.legend()
    # plt.savefig("../figs/ec1c_compare.png")

    figsetup(title="Testing for convergence of algorithm", xlab="x", ylab="u(x)",
             fname="ex1c_compare", show=showplots)

    return


def ex_d(showplots=False):

    # Require confirm because function takes a while to run
    confirm = raw_input("Do you want to run timing script? [y/n]:")

    if confirm.lower() != "y":
        print "Exiting program"
        sys.exit()

    general_times = []
    special_times = []

    List_of_N = [10, 25, 50, 75, 100, 250, 500, 750, 1e3, 2.5e3, 5e3, 7.5e3,
                 1e4, 1e5, 5e5, 1e6, 1e7]

    for N in List_of_N:
        N = int(N)               # Float input causes errors
        # Preparing input vectors before function call
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1      # Above diagonal

        print "\nStarting timing for N=%.1E" % N

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
             ylab="Time [s]", fname="ex1d_time", show=showplots)
    # Percentage difference in timings
    # Convert to arrays for easier manipulation
    general_times = np.array(general_times)
    special_times = np.array(special_times)

    percent_diff = np.abs(general_times - special_times) / \
        ((special_times + general_times) / 2.0) * 100

    plt.figure(figsize=(3.8, 3.8))
    plt.plot(log10N, percent_diff, "x--")
    figsetup(title="Percentage difference of execution time", xlab="$log_{10}N$",
             ylab="Percentage difference [%]", fname="ex1d_timediff", show=showplots)

    return


def ex_e(showplots=False):
    "Function calls to 'solve' question e. should combine with (d) for better efficiency"

    print "Starting (e)"

    List_of_N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
    errors = []

    for N in List_of_N:
        print "N=%.1E" % N
        N = int(N)
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1      # Above diagonal
        x, u = gauss_general(N, a, b, c)
        #x, u = gauss_specialized(N)
        #x, u = LU_benchmark(N)
        u2 = analyticSolution(x)
        errors.append(relError(u, u2))

    errors = np.array(errors)
    steps = 1 / np.array(List_of_N)

    plt.figure(figsize=(3.8, 3.8))
    plt.plot(np.log10(steps), np.log10(errors), "x--")
    figsetup(title="Error of algorithm for different step sizes", xlab="$log_{10}h$",
             ylab="$log_{10}\\epsilon_i$", fname="ex1e_err3", show=showplots)

    return


if __name__ == '__main__':
    #ex_c(showplots=True)
    ex_d(showplots=True)
    # ex_e(True)
