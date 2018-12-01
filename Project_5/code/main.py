# Python V2.7.15
# Contains all the calls to generate plots and other data for the report

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time
import os

relpath = "../data/"  # Relative path to folder were data files are stored.


def figsetup(title, xlab, ylab, fname, legend=True, show=False, tightlayout=True):
    """
    Sets up and saves figure for usage in report
    usage:
    plt.figure(figsize=[x, y])
    plot(...)
    ...
    plot(...)
    figsetup(title="filename", ...)
    """
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(fname)
    plt.title(title)
    if tightlayout is True:
        plt.tight_layout()
    if legend is True:
        plt.legend()
    plt.savefig("../figs/" + fname + ".pdf")
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def part_a():
    S = np.load(relpath + "S.npy")
    I = np.load(relpath + "I.npy")
    R = np.load(relpath + "R.npy")

    plt.plot(S)
    plt.plot(I)
    plt.plot(R)
    plt.show()
    return


def main():
    part_a()
    return


if __name__ == '__main__':
    main()
