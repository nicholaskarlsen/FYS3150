from __future__ import division
import matplotlib.pyplot as plt


def figsetup(title, xlab, ylab, fname, legend=True, show=False, tightlayout=True):
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
