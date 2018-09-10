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
    List_of_N = [1e1, 1e2, 1e3]
    plt.figure(figsize=(3.8, 3.8))
    for N in List_of_N:
        print "N=%.1E" % N
        N = int(N)  # Float input causes errors
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1       # Above diagonal

        x_gauss, u_gauss = general(N, a, b, c)

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
    """
    This function is exceedingly badly written as i had some issues when trying to 
    time my code in a previous itteration of the script. Since time is short, i opted
    to not spend any time to re-write this to a cleaner form once i got it all working.
    """

    print "Starting (d)"

    numcalls = 10

    List_of_N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    gen_times = []
    spe_times = []
    lu_times = []

    for N in List_of_N:
        print "N = %.1E" % N
        N = int(N)

        if N > 1e4:
            numcalls = 100
        else:
            numcalls = 10000
        exectime_gen = np.zeros(numcalls)

        print "Calling general..."

        for i in range(numcalls):
            a = np.ones(N) * -1
            b = np.ones(N) * 2
            c = np.ones(N) * -1
            t0 = time.time()
            x, u = general(N, a, b, c)
            t1 = time.time()
            exectime_gen[i] = t1 - t0
        gen_times.append(np.mean(exectime_gen))

        print "Calling specialized..."

        exectime_spe = np.zeros(numcalls)
        for i in range(numcalls):
            t0 = time.time()
            x, u = specialized(N)
            t1 = time.time()
            exectime_spe[i] = t1 - t0
        spe_times.append(np.mean(exectime_spe))

        if N <= 1e4:  # N > 1e4 -> run out of memory

            print "Calling LU..."
            exectime_lu = np.zeros(numcalls)
            for i in range(numcalls):
                t0 = time.time()
                x, u = specialized(N)
                t1 = time.time()
            exectime_lu[i] = t1 - t0
            lu_times.append(np.mean(exectime_lu))

    List_of_N = np.array(List_of_N)  # So i can operate on the entire lists
    gen_times = np.array(gen_times)
    spe_times = np.array(spe_times)
    lu_times = np.array(lu_times)

    plt.figure(figsize=(3.8, 3.8))
    plt.plot(np.log10(List_of_N), np.log10(gen_times), "x--", label="General algorithm")
    plt.plot(np.log10(List_of_N), np.log10(spe_times), "x--", label="specialized algorithm")
    plt.plot(np.log10(List_of_N[List_of_N<=1e4]), np.log10(lu_times), "x--", label="LU-Decomposition")

    figsetup(title="Timing algorithms", xlab="$log_{10}N$", ylab="$log_{10}t$",
             fname="ex1d_time", show=showplots)

    rdiff_time = np.abs(gen_times - spe_times) / ((np.abs(gen_times) + np.abs(spe_times)) / 2.0)
    plt.figure(figsize=(3.8, 3.8))
    plt.plot(np.log10(List_of_N), rdiff_time, "x--")
    figsetup(title="Relative difference of execution time", xlab="$log_{10}N$", ylab="Relative difference",
             fname="ex1d_timediff", show=showplots)

    return


def ex_e(showplots=False):
    "Function calls to 'solve' question e. should combine with (d) for better efficiency"

    print "Starting (e)"

    List_of_N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
    errors = []

    errors_LU = []

    for N in List_of_N:
        print "N=%.1E" % N
        N = int(N)
        a = np.ones(N) * -1      # Below diagonal
        b = np.ones(N) * 2       # Diagonal entries
        c = np.ones(N) * -1      # Above diagonal
        x, u = general(N, a, b, c)
        # x, u = specialized(N)
        # x, u = LU_benchmark(N)
        u2 = analyticSolution(x)
        errors.append(relError(u, u2))

        if N <= 1e4:
            xlu, ulu = LU_benchmark(N)
            u3 = analyticSolution(xlu)
            errors_LU.append(relError(ulu, u3))

    errors = np.array(errors)
    errors_LU = np.array(errors_LU)
    steps = 1 / np.array(List_of_N)
    print len(steps[steps<= 1/1e4]), len(errors_LU)
    plt.figure(figsize=(3.8, 3.8))
    plt.plot(np.log10(steps), np.log10(errors), "x--", label="algorithm")
    plt.plot(np.log10(steps[steps<= 1/1e4]), np.log10(errors_LU), "x--", label="LU")
    figsetup(title="Error of algorithm for different step sizes", xlab="$log_{10}h$",
             ylab="$log_{10}\\epsilon_i$", fname="ex1e_err", show=showplots)

    return


if __name__ == '__main__':
    # ex_c()
    #ex_d()
    ex_e()
