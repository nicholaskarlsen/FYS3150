import numpy as np
import matplotlib.pyplot as plt
from thomas_algorithm import *
from gausselim import *

# Initial Conditions & parameters
n = 1e3
h = 1.0 / float(n)
u0 = 0  # u(0) = 0
u1 = 0  # u(1) = 0
x = np.linspace(0, 1, h)


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
    plt.legend
    plt.savefig("../figs/" + fname + ".png", dpi=250)
    if show == False:
        plt.close()
    else:
        plt.show()
    return


def analyticSolution(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


if __name__ == '__main__':

    xvals = np.linspace(0, 1, 1e3)

    plt.plot(xvals, analyticSolution(xvals), label="Analytic Solution")
    for num in [10, 100, 1000]:
        plt.plot(gauss_general(num, f_func)[0], gauss_general(num, f_func)[
                 1], "x-", label="Tridiagonal Solution (n=%i)" % num)
    figsetup("test", "x", "y", "testfile")
