from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import metropolis as met
import time
import os


def figsetup(title, xlab, ylab, fname, legend=True, show=False, tightlayout=True):
    """
    Sets up and saves figure for usage in report
    usage:
    plt.figure(figsize=[x, y])
    plot(...)
    ...
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


def theoreticalVals(T):
    # Assumes J = 1
    # Expectation values
    E_theo = -8 * np.sinh(8 / T) / (3 + np.cosh(8 / T))
    E2_theo = (64 * np.cosh(8 / T)) / (np.cosh(8 / T) + 3)
    M_theo = 0
    M_abs_theo = (2 * np.exp(8 / T) + 4) / (np.cosh(8 / T) + 3)
    M2_abs_theo = (8 * np.exp(8 / T) + 8) / (np.cosh(8 / T) + 3)
    # Variance
    E_var = E2_theo - E_theo**2
    M_abs_var = M2_abs_theo - M_abs_theo**2
    # Derived variables
    C_v_theo = E_var / T**2
    chi_theo = M_abs_var / T
    return np.array([E_theo, M_theo, M_abs_theo, C_v_theo, chi_theo])


def prob_b():
    # Computes expectation values for different number of monte-carclo cycles
    # Many times to check consistency / error of results

    print "Running 2x2 test, may take a while"

    N = 2  # Size of lattice
    T = 1  # Temperature
    initState_o = met.generateState(N, 1)  # Ordered State
    initState_d = met.generateState(N)     # Dissordered State
    cycles = np.array([
        1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7
    ])
    n_samples = 100  # Number of times run algorithm per number of cycles
    # Initialize arrays to store results
    E = np.zeros([len(cycles), n_samples])
    M = np.zeros([len(cycles), n_samples])
    M_abs = np.zeros([len(cycles), n_samples])
    Cv = np.zeros([len(cycles), n_samples])
    chi = np.zeros([len(cycles), n_samples])

    # Setting up filenames for the results
    fnames = ["E_2x2bench.npy", "M_2x2bench.npy", "M_abs_2x2bench.npy",
              "cv_2x2bench.npy", "chi_2x2bench.npy"]
    for i in range(len(fnames)):
        fnames[i] = "conv_%i_" % n_samples + fnames[i]

    relpath = "../data/"

    # Check if file is saved
    if os.path.isfile(relpath + fnames[0]) is False:
        # if it is not, run algorithm (Which takes ~30-60 mins)
        for i in range(len(cycles)):              # For each cycles
            for j in range(n_samples):            # Run algorithm n_samples times & store results
                E[i, j], M[i, j], M_abs[i, j], Cv[i, j], chi[i, j] = met.montecarlo(spins=initState_o, T=1, trials=cycles[i])[-1]
                print "Running: No. Trials - %.1e, Sample No. - %i/%i, <E> = %.3f, <M> = %.3f " % (cycles[i], j, n_samples, E[i, j], M[i, j])
        # Save data for future use
        np.save(relpath + fnames[0], E)
        np.save(relpath + fnames[1], M)
        np.save(relpath + fnames[2], M_abs)
        np.save(relpath + fnames[3], Cv)
        np.save(relpath + fnames[4], chi)
    else:
        # Load data
        E = np.load(relpath + fnames[0])
        M = np.load(relpath + fnames[1])
        M_abs = np.load(relpath + fnames[2])
        Cv = np.load(relpath + fnames[3])
        chi = np.load(relpath + fnames[4])
        print "Found data, for n=%i samples" % len(E[0])

    theoVals = theoreticalVals(T) / N**2

    # Compute standard deviation
    E_std = [np.std(E[i]) for i in range(len(cycles))]
    M_std = [np.std(M[i]) for i in range(len(cycles))]
    M_abs_std = [np.std(M_abs[i]) for i in range(len(cycles))]
    Cv_std = [np.std(Cv[i]) for i in range(len(cycles))]
    chi_std = [np.std(chi[i]) for i in range(len(cycles))]

    # Generate plots
    width = 3
    height = 5

    # Plots for expectation values produced by algorithm for each parameter.
    # Energy
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[0], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), E[i], "x", color="blue")
    figsetup(title="Computed $\\left< E \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left< E \\right>/L^2$",
             fname="exb_convergencetest_E", legend=False)
    # Magnetic moment
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[1], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), M[i], "x", color="blue")
    figsetup(title="Computed $\\left< M \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left<M \\right>/L^2$",
             fname="exb_convergencetest_M", legend=False)
    # Magnetization
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[2], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), M_abs[i], "x", color="blue")
    figsetup(title="Computed $\\left< |M| \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left< |M| \\right>/L^2$",
             fname="exb_convergencetest_M_abs", legend=False)
    # Speciffic Heat Capacity
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[3], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), Cv[i], "x", color="blue")
    figsetup(title="Computed $C_v/L^2$", xlab="No. Monte Carlo cycles", ylab="$C_v/L^2$",
             fname="exb_convergencetest_Cv", legend=False)
    # Magnetic Suceptibility
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[4], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), chi[i], "x", color="blue")
    figsetup(title="Computed $\\chi/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\chi/L^2$",
             fname="exb_convergencetest_chi", legend=False)
    # Standard deviation of expectation values produced by algorithm
    # Energy
    plt.figure(figsize=[width, height])
    plt.loglog(cycles, E_std, "o--")
    figsetup(title="Std. of computed $\\left< E \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="Std. of Results",
             fname="exb_convergencetest_E_std", legend=False)
    # Magnetic moment
    plt.figure(figsize=[width, height])
    plt.loglog(cycles, M_std, "o--")
    figsetup(title="Std. of computed $\\left< M \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="Std. of Results",
             fname="exb_convergencetest_M_std", legend=False)
    # Magnetization
    plt.figure(figsize=[width, height])
    plt.loglog(cycles, M_abs_std, "o--")
    figsetup(title="Std. of computed $\\left< |M| \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="Std. of Results",
             fname="exb_convergencetest_M_abs_std", legend=False)
    # Speciffic heat capacity
    plt.figure(figsize=[width, height])
    plt.loglog(cycles, Cv_std, "o--")
    figsetup(title="Std. of computed $C_v/L^2$", xlab="No. Monte Carlo cycles", ylab="Std. of Results",
             fname="exb_convergencetest_Cv_std", legend=False)
    # Magnetic Suceptibility
    plt.figure(figsize=[width, height])
    plt.loglog(cycles, chi_std, "o--")
    figsetup(title="Std. of computed $\\chi/L^2$", xlab="No. Monte Carlo cycles", ylab="Std. of Results",
             fname="exb_convergencetest_chi_std", legend=False)
    # Box plot
    """
    plt.boxplot(E[:, -3:])
    plt.xticks(range(1, len(cycles[-3:]) + 1), ["%.e" % val for val in cycles[-3:]])
    figsetup(title="", xlab="No. Monte Carlo cycles", ylab="$\\left< E \\right>$",
             fname="exb_convergencetest_boxplot_E", legend=False)
    """
    print "2x2 test ----- Done"

    return


def prob_c():
    n_samples = 100  # number of samples per cycle
    T = [1, 2.4]

    sec1 = np.linspace(1, 100, 100)
    sec2 = np.linspace(int(1e2), int(1e3), 20)
    sec3 = np.linspace(int(1e3), int(1e4), 10)
    mc_trials = np.append(sec1, sec2, axis=0)
    mc_trials = np.append(mc_trials, sec3, axis=0)
    num_trials = len(mc_trials)

    E_mean_o = np.zeros([2, num_trials])
    E_mean_d = np.zeros([2, num_trials])
    M_mean_o = np.zeros([2, num_trials])
    M_mean_d = np.zeros([2, num_trials])

    # needs to store at most num_trials values in third nested array
    counter_o = np.zeros([2, num_trials])
    counter_d = np.zeros([2, num_trials])

    relpath = "../data/"
    fnames = [relpath + "prob_c_" + var + ".npy" for var in ["E_d_", "E_o_", "M_d_,", "M_o_", "c_d", "c_o"]]

    if os.path.isfile(fnames[0]) is False:
        # First ordered states
        initState_o = met.generateState(20, -1)  # Re-use same ordered state for all computations
        for i in range(len(T)):
            for j, trials in enumerate(mc_trials):  # For each number of monte carlo cycles
                print "Running trials:", j
                for k in range(n_samples):  # Compute mean 100 times and return mean result
                    # returns E, M, [E_mean, M_mean, M_abs_mean, Cv, chi], counter
                    tmp_o = met.montecarlo(spins=initState_o, T=T[i], trials=trials)
                    E_mean_o[i][j] += tmp_o[2][0]
                    M_mean_o[i][j] += tmp_o[2][2]
                    counter_o[i][j] += tmp_o[-1]
                    # Dissordered state, with new one each sample!
                    tmp_d = met.montecarlo(spins=met.generateState(20), T=T[i], trials=trials)
                    E_mean_d[i, j] += tmp_d[2][0]
                    M_mean_d[i, j] += tmp_d[2][2]
                    counter_d[i, j] += tmp_d[-1]
                # Normalize
                E_mean_o[i][j] /= n_samples
                M_mean_o[i][j] /= n_samples
                counter_o[i][j] /= n_samples
                E_mean_d[i][j] /= n_samples
                M_mean_d[i][j] /= n_samples
                counter_d[i][j] /= n_samples

        np.save(fnames[0], E_mean_d)
        np.save(fnames[1], E_mean_o)
        np.save(fnames[2], M_mean_d)
        np.save(fnames[3], M_mean_o)
        np.save(fnames[4], counter_d)
        np.save(fnames[5], counter_o)

    else:
        E_mean_d = np.load(fnames[0])
        E_mean_o = np.load(fnames[1])
        M_mean_d = np.load(fnames[2])
        M_mean_o = np.load(fnames[3])
        counter_d = np.load(fnames[4])
        counter_o = np.load(fnames[5])

    width = 4
    height = 2.5

    for i in range(2):
        plt.figure(figsize=[width, height])
        plt.plot(mc_trials, E_mean_o[i], "rx", label="Ordered", alpha=0.5)
        plt.plot(mc_trials, E_mean_d[i], "bx", label="Dissordered", alpha=0.5)
        figsetup(title="$k_B$T = %.1f" % T[i], xlab="No. Monte Carlo cycles", ylab="Mean computed $\\left<E\\right>/L^2$",
                 fname="exc_E_%i" % T[i], legend=True)

        plt.figure(figsize=[width, height])
        plt.plot(mc_trials, M_mean_o[i], "rx", label="Ordered", alpha=0.5)
        plt.plot(mc_trials, M_mean_d[i], "bx", label="Dissordered", alpha=0.5)
        figsetup(title="$k_B$T = %.1f" % T[i], xlab="No. Monte Carlo cycles", ylab="Mean computed $\\left<|M|\\right>/L^2$",
                 fname="exc_M_%i" % T[i], legend=True)

        plt.figure(figsize=[width, height])
        plt.plot(mc_trials, counter_o[i], "rx", label="Ordered", alpha=0.5)
        plt.plot(mc_trials, counter_d[i], "bx", label="Dissordered", alpha=0.5)
        figsetup(title="$k_B$T = %.1f" % T[i], xlab="No. Monte Carlo cycles", ylab="No. accepted configs.",
                 fname="exc_count_%i" % T[i], legend=True)

    return


def prob_d():
    """
    Finds P(E) by extracting the energies from algorithm and
    plots them in a normalized histogram rather than adding counters etc to the
    algorithm itself.
    """
    initState = met.generateState(20)  # Random 20x20 state
    trials = 1e7

    relpath = "../data/"
    fnames = [relpath + "prob_d_" + var + ".npy" for var in ["E", "M"]]

    E = np.zeros([2, int(trials)])
    M = np.zeros([2, int(trials)])
    T = [1, 2.4]

    if os.path.isfile(fnames[0]) is False:
        for i in range(len(T)):
            E[i], M[i] = met.montecarlo(spins=initState, T=T[i], trials=trials)[:2]
        np.save(fnames[0], E)
        np.save(fnames[1], M)
    else:
        E = np.load(fnames[0])
        M = np.load(fnames[1])

    eq = int(1e3)  # Equilisation period

    # Normalize the plot
    def weight(myarray): return np.ones_like(myarray) / float(len(myarray))
    for i in range(2):
        var = np.var(E[i, eq:])  # Variance
        plt.figure(figsize=[4, 4])
        plt.hist(E[i, eq:] / 20**2, weights=weight(E[i, eq:]), bins=50)
        figsetup(title="$k_B$T = %.1f" % T[i], xlab="$E$", ylab="No. Relative occurances",
                 fname="ex_d_histo_%i" % T[i], legend=False)
    return


def prob_e():
    # Calls code relevant to problem e
    latSize = [40, 60, 80, 100]

    """Wish to perform two sweeps, one rough which covers 6 pts in the range 2->2.3
    and a fine sweep, covering 10pts in the range 2.2 -> 2.3, where the phase transition
    Is expected to occur"""

    # Rough sweep
    N1 = 6
    T1 = np.linspace(2.0, 2.3, N1)
    # Fine sweep
    N2 = 10
    T2 = np.linspace(2.24, 2.29, N2)
    T = np.zeros(N1 + N2)
    T[:N1] = T1
    T[-N2:] = T2

    mc_cycles = 2e6

    # Set up arrays for storing results
    E = np.zeros([len(T), len(latSize)])
    M = np.zeros([len(T), len(latSize)])
    M_abs = np.zeros([len(T), len(latSize)])
    Cv = np.zeros([len(T), len(latSize)])
    chi = np.zeros([len(T), len(latSize)])

    # Setting up filenames for the results
    relpath = "../data/"
    pre = "prob_e_" + "%i_" % mc_cycles
    post = ".npy"
    fnames = [pre + var + post for var in ["E", "M", "M_abs", "cv", "chi"]]

    # Check if file is saved
    if os.path.isfile(relpath + fnames[0]) is False:
        # if it is not, run algorithm (Which takes ~30-60 mins)
        for i in xrange(len(T)):              # For each cycles
            print "Computing T=%.2f" % T[i]
            for j in xrange(len(latSize)):
                print "-- L= %i" % latSize[j]
                initState = met.generateState(latSize[j])
                E[i, j], M[i, j], M_abs[i, j], Cv[i, j], chi[i, j] = met.montecarlo(spins=initState, T=T[i], trials=mc_cycles)[-1]
        # Save data for future use
        np.save(relpath + fnames[0], E)
        np.save(relpath + fnames[1], M)
        np.save(relpath + fnames[2], M_abs)
        np.save(relpath + fnames[3], Cv)
        np.save(relpath + fnames[4], chi)
    else:
        # Load data
        E = np.load(relpath + fnames[0])
        M = np.load(relpath + fnames[1])
        M_abs = np.load(relpath + fnames[2])
        Cv = np.load(relpath + fnames[3])
        chi = np.load(relpath + fnames[4])
        print "Found data, for n=%i cycles" % mc_cycles

    # Sort arrays in ascending order of T
    permute = T.argsort()
    T = T[permute]
    E = E[permute]
    M = M[permute]
    M_abs = M_abs[permute]
    Cv = Cv[permute]
    chi = chi[permute]

    width = 9
    height = 2.75

    T_C = 2.269

    plt.figure(figsize=[width, height])
    plt.axvline(T_C, color="black", alpha=0.3)
    for i in range(len(latSize)):
        plt.plot(T, E[:, i], "x--", label="L=%i" % latSize[i])
    plt.xlim(2.15, 2.3)
    # Adds marker for T_C to xticks
    plt.xticks(list(plt.xticks()[0]) + [T_C], list(["%.2f" % s for s in plt.xticks()[0]] + ["$T_C$"]))
    figsetup(title="", ylab="$E/L^2$", xlab="$k_B$T", fname="ex_e_E")

    plt.figure(figsize=[width, height])
    plt.axvline(T_C, color="black", alpha=0.3)
    for i in range(len(latSize)):
        plt.plot(T, M_abs[:, i], "x--", label="L=%i" % latSize[i])
    plt.xlim(2.15, 2.3)
    plt.xticks(list(plt.xticks()[0]) + [T_C], list(["%.2f" % s for s in plt.xticks()[0]] + ["$T_C$"]))
    figsetup(title="", ylab="$|M|/L^2$", xlab="$k_B$T", fname="ex_e_M_abs")

    plt.figure(figsize=[5, 2.5])
    plt.axvline(T_C, color="black", alpha=0.3)
    for i in range(len(latSize)):
        plt.plot(T, M[:, i], "x--", label="L=%i" % latSize[i])
    plt.xlim(2.15, 2.3)
    plt.xticks(list(plt.xticks()[0]) + [T_C], list(["%.2f" % s for s in plt.xticks()[0]] + ["$T_C$"]))
    figsetup(title="", ylab="$M/L^2$", xlab="$k_B$T", fname="ex_e_M")

    plt.figure(figsize=[width, height])
    plt.axvline(T_C, color="black", alpha=0.3)
    for i in range(len(latSize)):
        plt.plot(T, Cv[:, i], "x--", label="L=%i" % latSize[i])
    plt.xlim(2.15, 2.3)
    plt.xticks(list(plt.xticks()[0]) + [T_C], list(["%.2f" % s for s in plt.xticks()[0]] + ["$T_C$"]))
    figsetup(title="", ylab="$C_v/L^2$", xlab="$k_B$T", fname="ex_e_Cv")

    plt.figure(figsize=[width, height])
    plt.axvline(T_C, color="black", alpha=0.3)
    for i in range(len(latSize)):
        plt.plot(T, chi[:, i], "x--", label="L=%i" % latSize[i])
    plt.xlim(2.15, 2.3)
    plt.xticks(list(plt.xticks()[0]) + [T_C], list(["%.2f" % s for s in plt.xticks()[0]] + ["$T_C$"]))
    figsetup(title="", ylab="$\\chi/L^2$", xlab="$k_B$T", fname="ex_e_chi")

    return


def main():
    # Uncomment the parts you wish to run
    # Note that the program will take a LONG time to run in the absence of
    # saved numpy arrays, which are not on github because of size limitations.

    print theoreticalVals(1) / 4
    prob_b()
    prob_c()
    prob_d()
    prob_e()
    return


if __name__ == '__main__':
    main()
