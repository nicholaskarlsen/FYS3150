# Written for Python Version 2.7.15
# By Nicholas Karlsen
# Based on script from lecture notes
from __future__ import division, print_function  # Nobody expects the integer division
import numpy as np
from numba import njit


def generateState(N, spin=False):
    # function wrapper which returns [N, N] array of random values in [-1, 1]
    # unless the spin parameter is set, in which case it will return a filled array
    if spin is False:
        array = np.random.choice([-1, 1], size=(N, N))
    else:
        array = np.full(shape=(N, N), fill_value=spin)

    return array


@njit   # JIT with enforced nopython=True
def montecarlo(spins, T, trials, returnCounter=False):
    """
    N: size of square lattice
    initState: Initial state of system, [N, N] array with entries +/-1
    T: Temperature of the system
    trials: Number of montecarlo cycles
    """

    # Finds correct indices accounting for periodic bounds
    # (NOTE: Has to be nested inside of montecarlo func for jit nopython mode)
    def periodic(i, limit, add):
        return (i + limit + add) % limit

    N = len(spins)                             # Ensure that N is int
    trials = int(trials)                       # ...           is int
    counter = 0                                # For counting accepted flips

    # Initialize variables
    E_mean = 0      # Expectation value of Energy
    E2_mean = 0     # ... Energy^2
    M_mean = 0      # ... Magnetization
    M2_mean = 0     # ... Magnetization^2
    M_abs_mean = 0  # ... Absolute Magnetization
    M2_abs_mean = 0  # ... Absolute Magnetization

    E = np.zeros(trials)
    M = np.zeros(trials)

    # Compute initial energy of the system
    for j in xrange(N):
        for i in xrange(N):
            E[0] -= spins[i, j] * (spins[periodic(i, N, -1), j] + spins[i, periodic(j, N, 1)])

    # And initial Magnetization
    M[0] = np.sum(spins)

    # Pre-compute possible change of energy
    w = np.zeros(17, dtype=np.float64)
    for dE in xrange(-8, 9, 4):  # include +8
        w[dE + 8] = np.exp(-dE / T)

    # Begin Monte-Carlo trials
    for i in xrange(1, trials + 1):
        # Perform metropolis algorithm
        # Set current energy to previous final energy before sweeping lattice
        E[i] = E[i - 1]
        M[i] = M[i - 1]
        for s in xrange(int(N**2)):     # Loop through N^2 randomly chosen spin-sites
            # Generate random positions within lattice
            x = np.random.randint(0, N)
            y = np.random.randint(0, N)
            # Compute change of energy
            dE = 2 * spins[x, y] *\
                (spins[periodic(x, N, -1), y] +
                 spins[periodic(x, N, 1), y] +
                 spins[x, periodic(y, N, -1)] +
                 spins[x, periodic(y, N, 1)])

            if np.random.random() <= w[dE + 8]:
                # Accept!
                spins[x, y] *= -1
                E[i] += dE
                M[i] += 2 * spins[x, y]
                counter += 1

        # Slightly faster than taking sum outside, 0.3s for trials=1e7
        E_mean += E[i]      # Expectation value of Energy
        E2_mean += E[i]**2     # ... Energy^2
        M_mean += M[i]      # ... Magnetization
        M2_mean += M[i]**2     # ... Magnetization^2
        M_abs_mean += abs(M[i])  # ... Absolute Magnetization
        M2_abs_mean += abs(M[i])**2  # ... Absolute Magnetization

    # Compute mean values
    E_mean /= trials
    E2_mean /= trials
    M_mean /= trials
    M2_mean /= trials
    M_abs_mean /= trials
    M2_abs_mean /= trials

    # Compute variance and normalize to per-point and temp
    E_variance = (E2_mean - E_mean**2) / float(N**2 * T**2)
    #M_variance = (M2_mean - M_mean**2) / float(N**2 * T**2)
    M_abs_variance = (M2_abs_mean - M_abs_mean**2) / float(N**2 * T**2)
    # Normalize returned averages to per-point
    E_mean /= float(N**2)       # Mean energy
    M_mean /= float(N**2)       # Mean magnetic moment
    M_abs_mean /= float(N**2)   # Mean magnetization
    # Compute derived values
    Cv = E_variance / T**2      # Speciffic Heat Capacitance
    chi = M_abs_variance / T        # Magnetic Suceptibility

    return E, M, [E_mean, M_mean, M_abs_mean, Cv, chi], counter


def main():
    initspins = generateState(2, 1)
    montecarlo(spins=initspins, T=1, trials=1e6)
    return


if __name__ == '__main__':
    main()
