# Python V2.7.15
# Contains all the calls to generate plots and other data for the report

from __future__ import division     # No body expects the integer division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import time
import os
import scipy.interpolate
from julia import Main as jcall  # Imports the main namespace of julia

import SIRS_ODE
import SIRS_MC

# Colours used in plots
S_colour = "blue"
I_colour = "red"
R_colour = "green"
Mean_colour = "orange"
ODE_colour = "black"
SIR_alpha = 0.05


def export_legend(legend, filename="../figs/prob_ab_legend.pdf", expand=[-5, -5, 5, 5]):
    # full credits to some stack overflow post for this & most of common_legends()
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


def common_legends():
    # Generates common legend for plots as separate figure.
    # First set up throw away figure with legend properties
    fig = plt.figure()
    ax = plt.gca()
    plt.plot([1, 2], [1, 2], linestyle="-", color=S_colour,
             label="S")
    plt.plot([1, 2], [1, 2], linestyle="-", color=I_colour,
             label="I")
    plt.plot([1, 2], [1, 2], linestyle="-", color=R_colour,
             label="R")
    plt.plot([1, 2], [1, 2], linestyle="--",
             color=Mean_colour, label="Mean MC")
    plt.plot([1, 2], [1, 2], linestyle="-", color=ODE_colour, label="ODE")
    # generate separate legend with properties from faux figure

    legend = plt.legend(ncol=5)  # Make the legend expand horizontally

    export_legend(legend)
    plt.close()

    fig = plt.figure()
    ax = plt.gca()
    # Multiply normal alpha by 10 to make it readable in legend.
    plt.plot([1, 2], [1, 2], linestyle="-", color=S_colour,
             label="S")
    plt.plot([1, 2], [1, 2], linestyle="-", color=I_colour,
             label="I")
    plt.plot([1, 2], [1, 2], linestyle="-", color=R_colour,
             label="R")
    # generate separate legend with properties from faux figure
    legend = plt.legend(ncol=3)  # Make the legend expand horizontally

    export_legend(legend, filename="../figs/prob_b_std_legend.pdf")
    plt.close()

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
    for N in [10000, 1000, 100, 10]:
        inst = SIRS_ODE.SIRS(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, N=N, tN=10)
        inst.solve(inst.sirs_basic)
        t, S, I, R = inst.get()
        h = 10 / float(N)
        plt.plot(t, S, label="h=%.2f" % h)
    plt.legend()
    plt.show()

    return


def part_a_b():

    print "Running a & b"

    """
    Produces plots etc relevant to question a and b.

    Note : All of the following was done before i wrote the julia MC solver,
    but since i saved .npy files after running code, theres no point in changing this
    to use julia version.
    """

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
        print b_i
        s, i, r = steadyState(a=4, b=b_i, c=0.5)  # Analytic steady state
        s *= S0 + I0 + R0  # Scale up to match population size
        i *= S0 + I0 + R0
        r *= S0 + I0 + R0

        # ~~~ ODE SOLUTIONS ~~~ #
        inst = SIRS_ODE.SIRS(S0=S0, I0=I0, R0=R0, a=a,
                             b=b_i, c=c, N=1000, tN=stop)
        inst.solve(inst.sirs_basic)
        t_ODE, S_ODE, I_ODE, R_ODE = inst.get()

        # Main Figure
        #fig, ax = plt.subplots(figsize=figdim)

        #ax.plot(t_ODE, S_ODE, color=S_colour)
        #ax.plot(t_ODE, I_ODE, color=I_colour)
        #ax.plot(t_ODE, R_ODE, color=R_colour)

        def plot_ab_settings():
            "Calls that are common to generating both figs"
            ax.text(x=1, y=350, s="b=%i" % b_i)
            ax.set_ylim(-50, S0 + I0 + R0 + 50)
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

        # plot_ab_settings()  # Apply settings to ODE fig
        #plt.savefig("../figs/prob_a_varb_%i.pdf" % b_i)
        # Save as .png as well because pdf reader didnt like vector plots for large no. trials
        #plt.savefig("../figs/prob_a_varb_%i.png" % b_i)
        # plt.close()

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
            t_MC, S_MC, I_MC, R_MC = SIRS_MC.main(
                S0=S0, I0=I0, R0=R0, a=a, b=b_i, c=c, stop_time=stop, trials=100)
            np.save(fn_t_curr, t_MC)
            np.save(fn_S_curr, S_MC)
            np.save(fn_I_curr, I_MC)
            np.save(fn_R_curr, R_MC)
        else:
            # Load if they exist
            t_MC = np.load(fn_t_curr)
            S_MC = np.load(fn_S_curr)
            I_MC = np.load(fn_I_curr)
            R_MC = np.load(fn_R_curr)

        # First plot all trials

        ax.plot(t_MC, S_MC, color=S_colour, alpha=SIR_alpha)
        ax.plot(t_MC, I_MC, color=I_colour, alpha=SIR_alpha)
        ax.plot(t_MC, R_MC, color=R_colour, alpha=SIR_alpha)

        # Then ODE solutions
        ax.plot(t_ODE, S_ODE, color=ODE_colour)
        ax.plot(t_ODE, I_ODE, color=ODE_colour)
        ax.plot(t_ODE, R_ODE, color=ODE_colour)

        # Then plot mean values
        ax.plot(t_MC, np.mean(S_MC, axis=1), label="S",
                color=Mean_colour, linestyle="--")
        ax.plot(t_MC, np.mean(I_MC, axis=1), label="I",
                color=Mean_colour, linestyle="--")
        ax.plot(t_MC, np.mean(R_MC, axis=1), label="R",
                color=Mean_colour, linestyle="--")

        plot_ab_settings()  # Apply same settings to MCMC fig
        plt.savefig("../figs/prob_b_varb_%i.pdf" % b_i)
        # NOTE : .pdf figs cause MAJOR stutter when viewing report. likely due to num. of pts.
        plt.savefig("../figs/prob_b_varb_%i.png" % b_i, dpi=600)
        plt.close()

    # Plot standard deviation for each b
    for b_i in list_of_bi:
        fn_S_curr = fn_S + "%i.npy" % b_i
        fn_I_curr = fn_I + "%i.npy" % b_i
        fn_R_curr = fn_R + "%i.npy" % b_i
        fn_t_curr = fn_t + "%i.npy" % b_i

        t_MC = np.load(fn_t_curr)
        S_MC = np.load(fn_S_curr)
        I_MC = np.load(fn_I_curr)
        R_MC = np.load(fn_R_curr)

        plt.figure(figsize=[5, 2])

        S_std = np.std(S_MC, axis=1)
        I_std = np.std(I_MC, axis=1)
        R_std = np.std(R_MC, axis=1)

        plt.plot(t_MC, S_std,
                 label="S", color=S_colour)
        plt.plot(t_MC, I_std,
                 label="I", color=I_colour)
        plt.plot(t_MC, R_std,
                 label="R", color=R_colour)

        plt.xlabel("Time")
        plt.ylabel("Standard Deviation")
        plt.title("b=%i" % b_i)
        plt.xlim(0, 30)
        plt.savefig("../figs/prob_b_std_%i.pdf" % b_i)
        plt.close()

    return


def part_c():

    print "Running c"

    jcall.include("SIRS_MC_vitdyn.jl")  # Import Julia program

    S0 = 300
    I0 = 100
    R0 = 0
    a = 4
    b = 1
    c = 0.5
    d = [1, 1, 1, 1]
    d_I = [0, 1, 1, 2]
    e = [1, 1, 1.2, 1.2]
    stop_time = [10, 60, 30, 7]
    trials = 100

    # Used to store number of steps performed in MC algorithm
    num_steps = np.zeros(trials)

    for i in range(4):
        plt.figure(figsize=[5, 2.5])

        # Compute ODE solution (Need t_ODE to compute interpolations)
        inst = SIRS_ODE.SIRS(
            S0=S0, I0=I0, R0=R0, N=1000, tN=stop_time[i], a=a, b=b, c=c, d=d[i], d_I=d_I[i], e=e[i]
        )
        inst.solve(inst.sirs_vitdyn)
        t_ODE, S_ODE, I_ODE, R_ODE = inst.get()

        S_mean = np.zeros(len(t_ODE))
        I_mean = np.zeros(len(t_ODE))
        R_mean = np.zeros(len(t_ODE))

        # Compute MC solutions
        for j in range(trials):
            command = "SIRS_vitdyn(S0=%i, I0=%i, R0=%i, a=%.2f, b=%.2f, c=%.2f, d=%.2f, \
                d_I=%.2f, e=%.2f, stop_time=%i)" % (S0, I0, R0, a, b, c, d[i], d_I[i], e[i], stop_time[i])
            t, S, I, R = jcall.eval(command)

            """ Note: This way of using PyJulia may not be the most efficient/proper way?
                The package lacks any meaningfull documentation by the looks of it, so
                this way of doing things will have to do. """

            plt.plot(t, S, color=S_colour, alpha=SIR_alpha)
            plt.plot(t, I, color=I_colour, alpha=SIR_alpha)
            plt.plot(t, R, color=R_colour, alpha=SIR_alpha)

            interpolated_S = scipy.interpolate.interp1d(t, S, kind="linear")
            interpolated_I = scipy.interpolate.interp1d(t, I, kind="linear")
            interpolated_R = scipy.interpolate.interp1d(t, R, kind="linear")

            t_ODE_valid = t_ODE[t_ODE < t[-1]]

            # Get interpolations at same pts as t_ODE is solved
            S_mean[:len(t_ODE_valid)] += interpolated_S(t_ODE_valid)
            I_mean[:len(t_ODE_valid)] += interpolated_I(t_ODE_valid)
            R_mean[:len(t_ODE_valid)] += interpolated_R(t_ODE_valid)

            if i == 1:
                # Sufficient to look at length of data sets for only one of the cycles
                num_steps[j] = len(t)

            # TODO : Add interpolation at points t_ODE using scipy.interpolate.inerp1d()
        S_mean /= float(trials)
        I_mean /= float(trials)
        R_mean /= float(trials)

        # Add ODE solution on top of MC solutions
        plt.plot(t_ODE, S_ODE, color=ODE_colour)
        plt.plot(t_ODE, I_ODE, color=ODE_colour)
        plt.plot(t_ODE, R_ODE, color=ODE_colour)

        # Dotted mean values on top of ODE
        plt.plot(t_ODE, S_mean, color=Mean_colour, linestyle="--")
        plt.plot(t_ODE, I_mean, color=Mean_colour, linestyle="--")
        plt.plot(t_ODE, R_mean, color=Mean_colour, linestyle="--")

        plt.xlim(0, stop_time[i])
        plt.ylim(-10, 450)
        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.title("a=%.2f, b=%.2f, c=%.2f, d=%.2f, $d_I$=%.2f, e=%.2f" %
                  (a, b, c, d[i], d_I[i], e[i]))
        plt.savefig("../figs/prob_c_fig_%i.pdf" % i)
        plt.savefig("../figs/prob_c_fig_%i.png" % i)
        plt.close()

    # Plot number of data points for cycles in one of the runs
    plt.figure(figsize=[5, 5])
    plt.hist(num_steps, bins=40)
    plt.xlabel("Number of steps in trial")
    plt.title("No. steps in 100 trials")
    plt.savefig("../figs/num_steps_c.png")
    plt.close()

    return


def part_d():

    print "Running d"

    jcall.include("SIRS_MC_svar.jl")  # Import Julia program

    # Keep following parameters constant
    S0 = 300
    I0 = 100
    R0 = 0
    a0 = 4
    b = 1
    c = 0.5
    # But change the following
    percent_diff = [1, 1, .5, .5]  # maximum percent deviation from a0
    # Amplitude of deviation from a0
    Amp = [diff * a0 for diff in percent_diff]
    omega = [4, 0.25, 2, 0.25]  # Frequency of deviation from a0
    stop = [10, 70, 10, 70]  # Simulation time

    trials = 100  # No. times to run MC simulation

    for i in range(len(percent_diff)):
        plt.figure(figsize=[5, 2.5])

        # Get MC data from julia program (SIRS_MC.jl)
        command = "SIRS_svar(S0=%i, I0=%i, R0=%i, a0=%.2f, A=%.2f, omega=%.2f, b=%.2f, c=%.2f, stop_time=%i, trials=%i)"\
            % (S0, I0, R0, a0, Amp[i], omega[i], b, c, stop[i], trials)
        t, S, I, R = jcall.eval(command)

        for j in range(trials):
            plt.plot(t, S[j], color=S_colour, alpha=SIR_alpha)
            plt.plot(t, I[j], color=I_colour, alpha=SIR_alpha)
            plt.plot(t, R[j], color=R_colour, alpha=SIR_alpha)

        # Plot ODE solution on top
        inst = SIRS_ODE.SIRS(
            S0=S0, I0=I0, R0=R0, N=int(1e3), tN=stop[i], a=a0, b=b, c=c, Amplitude=Amp[i],
            omega=omega[i])
        inst.solve(inst.sirs_svar)
        t_ODE, S_ODE, I_ODE, R_ODE = inst.get()
        plt.plot(t_ODE, S_ODE, color="Black")
        plt.plot(t_ODE, I_ODE, color="Black")
        plt.plot(t_ODE, R_ODE, color="Black")
        # Mean MC on top of ODE
        plt.plot(t, np.mean(S, axis=0), linestyle="--", color=Mean_colour)
        plt.plot(t, np.mean(I, axis=0), linestyle="--", color=Mean_colour)
        plt.plot(t, np.mean(R, axis=0), linestyle="--", color=Mean_colour)

        plt.xlim(0, stop[i])
        plt.ylim(-10, 410)
        plt.title("A=%.2f, $\\omega$=%.2f" % (Amp[i], omega[i]))
        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.savefig("../figs/prob_d_fig_%i.pdf" % i)
        plt.savefig("../figs/prob_d_fig_%i.png" % i)
        plt.close()

    return


def part_e():

    print "Running e"

    jcall.include("SIRS_MC_vax.jl")  # Include Julia program

    # Keep following parameters constant
    S0 = 300
    I0 = 100
    R0 = 0
    a = 4
    b = 1
    c = 0.5
    # But change the following
    f = [50, 100, 200, 300]
    stop = [20, 20, 20, 20]  # Simulation time
    trials = 100  # No. times to run MC simulation

    for i in range(len(f)):
        plt.figure(figsize=[5, 2.5])

        command = "SIRS_vax(S0=%i, I0=%i, R0=%i, a=%.2f, b=%.2f, c=%.2f, f=%.2f, stop_time=%i, trials=%i)"\
            % (S0, I0, R0, a, b, c, f[i], stop[i], trials)
        t_MC, S_MC, I_MC, R_MC = jcall.eval(command)
        for j in range(trials):
            plt.plot(t_MC, S_MC[j], color=S_colour, alpha=SIR_alpha)
            plt.plot(t_MC, I_MC[j], color=I_colour, alpha=SIR_alpha)
            plt.plot(t_MC, R_MC[j], color=R_colour, alpha=SIR_alpha)
        # Plot ODE solution on top
        inst = SIRS_ODE.SIRS(
            S0=S0, I0=I0, R0=R0, N=int(1e5), tN=stop[i], a=a, b=b, c=c, f=f[i]
        )
        inst.solve(inst.sirs_vax)
        t_ODE, S_ODE, I_ODE, R_ODE = inst.get()
        plt.plot(t_ODE, S_ODE, color=ODE_colour)
        plt.plot(t_ODE, I_ODE, color=ODE_colour)
        plt.plot(t_ODE, R_ODE, color=ODE_colour)
        # And finally MC mean on top of that.
        plt.plot(t_MC, np.mean(S_MC, axis=0),
                 linestyle="--", color=Mean_colour)
        plt.plot(t_MC, np.mean(I_MC, axis=0),
                 linestyle="--", color=Mean_colour)
        plt.plot(t_MC, np.mean(R_MC, axis=0),
                 linestyle="--", color=Mean_colour)

        plt.xlim(0, stop[i])
        plt.title("f=%.2f" % f[i])
        plt.xlabel("Time")
        plt.ylabel("No. People")
        plt.savefig("../figs/prob_e_fig_%i.pdf" % i)
        plt.savefig("../figs/prob_e_fig_%i.png" % i)
        plt.close()

        # NOTE: System poorly modeled by ODE because the "f" term is not
        # Population conservative. Will continue to "transfer" from
        # I -> S even if I >= 0

    return


def main():
    # convergence_check()
    common_legends()
    # part_a_b()
    part_c()
    # part_d()
    # part_e()
    return


if __name__ == '__main__':
    main()
