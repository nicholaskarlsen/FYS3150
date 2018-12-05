# Python V2.7.15
# Contains all the calls to generate plots and other data for the report

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import time
import os

import SIRS_ODE
import SIRS_MCMC


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
        plt.legend(loc="best")
    plt.savefig("../figs/" + fname + ".pdf")
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def steadyState(a, b, c):
    # Steady state for the "sirs_basic" method in SIRS_ODE
    s = b / a
    i = (1 - b / a) / (1 + b / c)
    r = b / c * (1 - b / a) / (1 + b / c)
    return s, i, r


def convergence_check():
    """
    Check how well the algorithm conserves population for different step sizes
    """
    for N in [1000, 100, 10]:
        inst = SIRS_ODE.SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=N, tN=10)
        inst.solve(inst.sirs_basic)
        t, S, I, R = inst.get()
        h = 10 / float(N)
        plt.plot(t, S, label="h=%.2f" % h)
    plt.legend()
    plt.show()

    return


def part_a_b():
    S0 = 300
    I0 = 100
    R0 = 0
    a = 4
    c = 0.5
    stop = 15

    figdim = [4, 2.3]
    for b_i in [1, 2, 3, 4]:
        s, i, r = steadyState(a=4, b=b_i, c=0.5)  # Analytic steady state
        s *= S0 + I0 + R0  # Scale up to match population size
        i *= S0 + I0 + R0
        r *= S0 + I0 + R0

        # ~~~ ODE SOLUTIONS ~~~ #
        inst = SIRS_ODE.SIRS(S0=S0, I0=I0, R0=R0, a=a, b=b_i, c=c, N=1000, tN=stop)
        inst.solve(inst.sirs_basic)
        t, S, I, R = inst.get()

        fig, ax = plt.subplots(figsize=figdim)
        # Main Figure
        def tmp():
            ax.plot(t, S, label="S", color="blue")
            ax.plot(t, I, label="I", color="red")
            ax.plot(t, R, label="R", color="green")
            ax.text(x=1, y=350, s="b=%i" % b_i)
            ax.set_ylim(0, S0 + I0 + R0 + 10)
            ax.set_xlim(0, stop)
            ax.set_xlabel("Time")
            ax.set_ylabel("No. People")
            ax.legend(loc="upper right")
            # Right y-axis s* i* r* ticks
            ax_ticks = ax.twinx()
            ax_ticks.figure.canvas.draw()
            y1, y2 = ax.get_ylim()
            ax_ticks.set_ylim(y1, y2)
            ax_ticks.set_yticks([s, i, r])
            ax_ticks.set_yticklabels(["$s^*$", "$i^*$", "$r^*$"])
            # Save
            plt.tight_layout()
        tmp()
        plt.savefig("../figs/prob_a_varb_%i.pdf" % b_i)
        plt.close()

        # ~~~ MCMC SOLUTION ~~~ #
        t, S, I, R = SIRS_MCMC.main(S0=S0, I0=I0, R0=R0, a=a, b=b_i, c=c, stop_time=stop)
        """
        plt.plot(t, S, label="S", color="blue")
        plt.plot(t, I, label="I", color="red")
        plt.plot(t, R, label="R", color="green")
        plt.text(x=1, y=350, s="b=%i" % b_i)
        plt.ylim(0, S0 + I0 + R0 + 10)
        plt.xlim(0, stop)
        figsetup("MCMC", "Time", "No. People", "prob_b_varb_%i" % b_i, legend=False)
        """
        fig, ax = plt.subplots(figsize=figdim)

        tmp()
        plt.savefig("../figs/prob_b_varb_%i.pdf" % b_i)
        plt.close()


    return


def main():
    # convergence_check()
    part_a_b()
    return


if __name__ == '__main__':
    main()
