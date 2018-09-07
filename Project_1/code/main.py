from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
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
    plt.figure(figsize=(3.8, 3.8))
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


if __name__ == '__main__':
    steps = []
    maxerrors = []

    list_of_N = np.linspace(1e2, 1e8, 50)
    for i in list_of_N:
        print "Computing for N=%i ..." %i
        i = int(i)
        x, u = gauss_general(i)
        h = 1 / i

        maxerrors.append(relError(u, analyticSolution(x)))
        steps.append(h)



    plt.plot(np.log10(steps), np.log10(maxerrors), "x--")
    plt.xlabel("Step size $log_{10}(\\Delta x)$")
    plt.ylabel("Max error $log_{10}\\epsilon(\\Delta x)$")
    plt.savefig("../figs/" + "error_general" + ".png")
    plt.show()