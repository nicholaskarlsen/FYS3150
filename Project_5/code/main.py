# Python V2.7.15
# Contains all the calls to generate plots and other data for the report

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import time
import os
from julia import Main as jcall  # Imports the main namespace of julia

import SIRS_ODE
import SIRS_MCMC

jcall.include("SIRS_MCMC.jl")  # Imports the file

# Colours used in plots
S_colour = "blue"
I_colour = "red"
R_colour = "green"
Mean_colour = "black"


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
    stop = 30

    fn_S = "../data/prob_b_S_"
    fn_I = "../data/prob_I_I_"
    fn_R = "../data/prob_b_R_"
    fn_t = "../data/prob_b_t_"

    figdim = [5, 2.5]

    list_of_bi = np.array([1, 2, 3, 4])

    variance = np.zeros(len(list_of_bi))  # Store variance of result
    # Alternatively look at difference between mean and  s*, i*, ... ?

    for b_i in list_of_bi:
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

        ax.plot(t, S, color=S_colour)
        ax.plot(t, I, color=I_colour)
        ax.plot(t, R, color=R_colour)

        def plot_ab_settings():
            "Calls that are common to generating both figs"
            ax.text(x=1, y=350, s="b=%i" % b_i)
            ax.set_ylim(0, S0 + I0 + R0 + 10)
            ax.set_xlim(0, stop)
            ax.set_xlabel("Time")
            ax.set_ylabel("No. People")
            # ax.legend(loc="upper right")
            # Right y-axis s* i* r* ticks
            ax_ticks = ax.twinx()
            ax_ticks.figure.canvas.draw()
            y1, y2 = ax.get_ylim()
            ax_ticks.set_ylim(y1, y2)
            ax_ticks.set_yticks([s, i, r])
            ax_ticks.set_yticklabels(["$s^*$", "$i^*$", "$r^*$"])
            # Save
            plt.tight_layout()

        plot_ab_settings()  # Apply settings to ODE fig
        plt.savefig("../figs/prob_a_varb_%i.pdf" % b_i)
        # Save as .png as well because pdf reader didnt like vector plots for large no. trials
        plt.savefig("../figs/prob_a_varb_%i.png" % b_i)
        plt.close()

        # ~~~ MCMC SOLUTION ~~~ #
        N_samples = 100

        # Filenames for current b
        fn_S_curr = fn_S + "%i.npy" % b_i
        fn_I_curr = fn_I + "%i.npy" % b_i
        fn_R_curr = fn_R + "%i.npy" % b_i
        fn_t_curr = fn_t + "%i.npy" % b_i

        # Initialize figure
        fig, ax = plt.subplots(figsize=figdim)

        # If one of the files dont exist, all of them probably dont
        if os.path.isfile(fn_S_curr) is False:
            # If they dont exist, compute & save
            t, S, I, R = SIRS_MCMC.main(S0=S0, I0=I0, R0=R0, a=a, b=b_i, c=c, stop_time=stop, trials=100)
            np.save(fn_t_curr, t)
            np.save(fn_S_curr, S)
            np.save(fn_I_curr, I)
            np.save(fn_R_curr, R)
        else:
            # Load if they exist
            t = np.load(fn_t_curr)
            S = np.load(fn_S_curr)
            I = np.load(fn_I_curr)
            R = np.load(fn_R_curr)

        # First plot all trials
        plt.plot(t, S, color=S_colour, alpha=.1)
        plt.plot(t, I, color=I_colour, alpha=.1)
        plt.plot(t, R, color=R_colour, alpha=.1)

        # Then plot mean values
        plt.plot(t, np.mean(S, axis=1), label="S", color="black")
        plt.plot(t, np.mean(I, axis=1), label="I", color="black")
        plt.plot(t, np.mean(R, axis=1), label="R", color="black")

        plot_ab_settings()  # Apply same settings to MCMC fig
        plt.savefig("../figs/prob_b_varb_%i.pdf" % b_i)
        plt.savefig("../figs/prob_b_varb_%i.png" % b_i)
        plt.close()

    # --- Generate separate figure for common legend --- #
    # First set up throw away figure with legend properties
    fig = plt.figure()
    ax = plt.gca()
    plt.plot([1, 2], [1, 2], linestyle="-", color=S_colour, label="S")
    plt.plot([1, 2], [1, 2], linestyle="-", color=I_colour, label="I")
    plt.plot([1, 2], [1, 2], linestyle="-", color=R_colour, label="R")
    plt.plot([1, 2], [1, 2], linestyle="-", color=Mean_colour, label="Mean")
    # generate separate legend with properties from faux figure

    legend = plt.legend(ncol=4)

    def export_legend(legend, filename="../figs/prob_ab_legend.pdf", expand=[-5, -5, 5, 5]):
        fig = legend.figure
        fig.canvas.draw()
        bbox = legend.get_window_extent()
        bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(filename, dpi="figure", bbox_inches=bbox)

    export_legend(legend)
    plt.close()

    # plot variace against dt

    return


"""
def part_c():

    stop = 20

    no_sets = 4
    no_trials = 100
    fn_set = range(1, no_sets + 1)
    fn_trials = range(1, no_trials + 1)

    for i, dset in enumerate(fn_set):
        plt.figure(figsize=[4, 2.3])
        for j, trial in enumerate(fn_trials):
            S_mc = np.load("../data/prob_c/c_%i/S_%i.npy" % (dset, trial))
            I_mc = np.load("../data/prob_c/c_%i/I_%i.npy" % (dset, trial))
            R_mc = np.load("../data/prob_c/c_%i/R_%i.npy" % (dset, trial))
            t_mc = np.load("../data/prob_c/c_%i/t_%i.npy" % (dset, trial))

            plt.plot(t_mc, S_mc, color="blue", alpha=.1)
            plt.plot(t_mc, I_mc, color="red", alpha=.1)
            plt.plot(t_mc, R_mc, color="green", alpha=.1)

        inst = SIRS_ODE.SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, e=1, d=1, d_I=i, N=100, tN=20)
        inst.solve(inst.sirs_basic)
        t_ode, S_ode, I_ode, R_ode = inst.get()

        plt.plot(t_ode, S_ode, color="black")
        plt.plot(t_ode, I_ode, color="black")
        plt.plot(t_ode, R_ode, color="black")
        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.xlim(0, stop)
        plt.tight_layout()
        plt.savefig("../figs/prob_c_fig_%i.pdf" % i)
        plt.savefig("../figs/prob_c_fig_%i.png" % i)
        plt.close()

    return
"""


def part_c():
    S0 = 300
    I0 = 100
    R0 = 0
    a = 4
    b = 1
    c = 0.5
    d = [1, 1, 1, 1]
    d_I = [1, 2, 0, 1]
    e = [1, 1, 1, 2]
    stop_time = 20
    trials = 100

    for i in range(4):
        plt.figure(figsize=[5, 2.5])
        for j in range(trials):
            command = "SIRS_vitdyn(S0=%i, I0=%i, R0=%i, a=%i, b=%i, c=%i, d=%i, d_I=%i,\
                e=%i, stop_time=%i)" % (S0, I0, R0, a, b, c, d[i], d_I[i], e[i], stop_time)
            t, S, I, R = jcall.eval(command)

            plt.plot(t, S, color=S_colour, alpha=0.1)
            plt.plot(t, I, color=I_colour, alpha=0.1)
            plt.plot(t, R, color=R_colour, alpha=0.1)
        # Add ODE solution
        inst = SIRS_ODE.SIRS(
            S0=S0, I0=I0, R0=R0, N=1000, tN=stop_time, a=a, b=b, c=c, d=d[i], d_I=d_I[i], e=e[i]
        )
        inst.solve(inst.sirs_basic)
        t, S, I, R = inst.get()
        plt.plot(t, S, color="Black")
        plt.plot(t, I, color="Black")
        plt.plot(t, R, color="Black")
        plt.xlim(0, stop_time)
        plt.ylim(0, 400)
        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.title("a=%i, b=%i, c=%.1f, d=%i, $d_I$=%i, e=%i" % (a, b, c, d[i], d_I[i], e[i]))
        plt.savefig("../figs/prob_c_fig_%i.pdf" % i)
        plt.savefig("../figs/prob_c_fig_%i.png" % i)
        plt.close()

    return


def part_d():
    S0 = 300
    I0 = 100
    R0 = 0
    a0 = 4
    b = 1
    c = 0.5

    percent_diff = [0.25, 0.25, 1, 1]
    Amp = [diff * a0 for diff in percent_diff]
    omega = [2, 0.2, 2, 0.5]
    stop_time = 20
    trials = 100

    for i in range(4):
        plt.figure(figsize=[5, 2.5])
        for j in range(trials):
            command = "SIRS_svar(S0=%i, I0=%i, R0=%i, a0=%i, A=%.2f, omega=%.2f, b=%i, c=%.2f, stop_time=%i)"\
                % (S0, I0, R0, a0, Amp[i], omega[i], b, c, stop_time)
            t, S, I, R = jcall.eval(command)

            plt.plot(t, S, color=S_colour, alpha=0.1)
            plt.plot(t, I, color=I_colour, alpha=0.1)
            plt.plot(t, R, color=R_colour, alpha=0.1)
        # Add ODE solution
        inst = SIRS_ODE.SIRS(
            S0=S0, I0=I0, R0=R0, N=int(1e3), tN=stop_time, a=a0, b=b, c=c, Amplitude=Amp[i],
            omega=omega[i]
        )
        inst.solve(inst.sirs_svar)
        t, S, I, R = inst.get()
        plt.plot(t, S, color="Black")
        plt.plot(t, I, color="Black")
        plt.plot(t, R, color="Black")
        plt.xlim(0, stop_time)
        plt.ylim(0, 400)
        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.savefig("../figs/prob_d_fig_%i.pdf" % i)
        plt.savefig("../figs/prob_d_fig_%i.png" % i)
        plt.close()

    return


def main():
    # convergence_check()
    # part_a_b()
    # part_b()
    #part_c()
    part_d()
    # part_e()
    return


if __name__ == '__main__':
    main()
