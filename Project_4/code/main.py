from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import metropolis as met
from scipy.constants import k
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
    E2_theo = - (64 * np.cosh(8 / T)) / (np.cosh(8 / T) + 3)
    M_theo = 0
    M_abs_theo = (2 * np.exp(8 / T) + 4) / (np.cosh(8 / T) + 3)
    M2_abs_theo = (4 * np.exp(8 / T)) / (np.cosh(8 / T) + 3)
    # Variance
    E_var = E2_theo - E_theo**2
    M_abs_var = M2_abs_theo - M_abs_theo**2
    # Derived variables
    C_v_theo = E_var / T**2
    chi_theo = M_abs_var / T

    return np.array([E_theo, M_theo, M_abs_theo, C_v_theo, chi_theo])


def prob_2():

    L = 25
    n_temps = 10
    temps = np.linspace(2.2, 2.4, n_temps)
    initstate = met.generateState(L, 1)
    E = np.zeros(n_temps)
    E2 = np.zeros(n_temps)
    for i in range(n_temps):
        E[i], E2[i] = met.montecarlo(spins=initstate, T=temps[i], trials=1e6)

    plt.plot(temps, E / L**2)
    plt.show()

    return


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
                E[i, j], M[i, j], M_abs[i, j], Cv[i, j], chi[i, j] = met.montecarlo(spins=initState_o, T=1, trials=cycles[i])
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

    # Mean difference from expectation value

    ##################
    # Generate plots #
    ##################
    # First plots for expectation values produced by algorithm for each parameter.
    # Energy

    width = 3
    height = 5

    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[0], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), E[i], "x", color="blue")
    figsetup(title="Computed $\\left< E \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left< E \\right>$",
             fname="exb_convergencetest_E", legend=False)
    # Magnetic moment
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[1], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), M[i], "x", color="blue")
    figsetup(title="Computed $\\left< M \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left<M \\right>$",
             fname="exb_convergencetest_M", legend=False)
    # Magnetization
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[2], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), M_abs[i], "x", color="blue")
    figsetup(title="Computed $\\left< |M| \\right>/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left< |M| \\right>$",
             fname="exb_convergencetest_M_abs", legend=False)
    # Speciffic Heat Capacity
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[3], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), Cv[i], "x", color="blue")
    figsetup(title="Computed $C_v/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left< |M| \\right>$",
             fname="exb_convergencetest_Cv", legend=False)
    # Magnetic Suceptibility
    plt.figure(figsize=[width, height])
    plt.axhline(theoVals[4], color="red", label="Theoretical Value", linestyle="--")
    for i in range(len(cycles)):
        plt.semilogx(cycles[i] * np.ones(n_samples), chi[i], "x", color="blue")
    figsetup(title="Computed $\\chi/L^2$", xlab="No. Monte Carlo cycles", ylab="$\\left< |M| \\right>$",
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

    initState = met.generateState(20)  # Random 20x20 state

    return


def prob_e():

    latSize = [40, 60, 80, 100]
    # Rough sweep
    N1 = 6
    T1 = np.linspace(2.0, 2.3, N1)
    # Fine sweep
    N2 = 10
    T2 = np.linspace(2.24, 2.29, N2)
    T = np.zeros(N1 + N2)
    T[:N1] = T1
    T[-N2:] = T2
    print T

    mc_cycles = 1e6

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
                E[i, j], M[i, j], M_abs[i, j], Cv[i, j], chi[i, j] = met.montecarlo(spins=initState, T=T[i], trials=mc_cycles)
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

    width = 4
    height = 4

    plt.figure(figsize=[width, height])

    for i in range(len(latSize)):
        plt.plot(T, E[:, i], "x--", label="L=%i" % latSize[i])
    plt.legend()
    plt.show()

    for i in range(len(latSize)):
        plt.plot(T, M_abs[:, i], "x--", label="L=%i" % latSize[i])
    plt.legend()
    plt.show()

    for i in range(len(latSize)):
        plt.plot(T, Cv[:, i], "x--", label="L=%i" % latSize[i])
    plt.legend()
    plt.show()
    for i in range(len(latSize)):
        plt.plot(T, chi[:, i], "x--", label="L=%i" % latSize[i])
    plt.legend()
    plt.show()

    return


def main():
    #prob_b()
    prob_e()
    return


if __name__ == '__main__':
    main()
