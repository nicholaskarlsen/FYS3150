# Written for Python Version 2.7.15
# By Nicholas Karlsen
# Based on script from lecture notes
from __future__ import division, print_function  # Nobody expects the integer division
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt

def generateState(N, spin=False):
    # function wrapper which returns [N, N] array of random values in [-1, 1]
    # unless the spin parameter is set, in which case it will return a filled array
    if spin is False:
        array = np.random.choice([-1, 1], size=(N, N))
    else:
        array = np.full(shape=(N, N), fill_value=spin)

    return array


@njit
def montecarlo(spins, T, trials):
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

    # Initialize variables
    E = 0           # Energy
    E_mean = 0      # Expectation value of Energy
    E2_mean = 0     # ... Energy^2
    M_mean = 0      # ... Magnetization
    M2_mean = 0     # ... Magnetization^2
    M_abs_mean = 0  # ... Absolute Magnetization

    # Compute initial energy of the system
    for j in xrange(N):
        for i in xrange(N):
            E -= spins[i, j] * (spins[periodic(i, N, -1), j] + spins[i, periodic(j, N, 1)])

    # And initial Magnetization
    M = np.sum(spins)

    # Pre-compute possible change of energy
    w = np.zeros(17, dtype=np.float64)
    for dE in xrange(-8, 9, 4):  # include +8
        w[dE + 8] = np.exp(-dE / T)

    # Begin Monte-Carlo trials
    for i in xrange(trials):
        # Perform metropolis algorithm
        for s in xrange(int(N**2)):     # Loop through N^2 randomly chosen spin-sites
            x = np.random.randint(0, N)
            y = np.random.randint(0, N)

            dE = 2 * spins[x, y] *\
                (spins[periodic(x, N, -1), y] +
                 spins[periodic(x, N, 1), y] +
                 spins[x, periodic(y, N, -1)] +
                 spins[x, periodic(y, N, 1)])

            if np.random.random() <= w[dE + 8]:
                # Accept!
                spins[x, y] *= -1
                E += dE
                M += 2 * spins[x, y]

        E_mean += E
        E2_mean += E**2
        M_mean += M
        M2_mean += M**2
        M_abs_mean += abs(M)

    # Normalize
    E_mean /= float(trials)         # Mean energy
    E2_mean /= float(trials)        # ... energy^2
    M_mean /= float(trials)         # ... Magnetic moment
    M2_mean /= float(trials)        # ... Magnetic moment^2
    M_abs_mean /= float(trials)     # ... Magnetization
    # Calculate variance and normalize to per-point and temp
    E_variance = (E2_mean - E_mean**2) / float(N**2 * T**2)
    M_variance = (E2_mean - E_mean**2) / float(N**2 * T**2)
    # Normalize returned averages to per-point
    E_mean /= float(N**2)       # Mean energy
    M_mean /= float(N**2)       # Mean magnetic moment
    M_abs_mean /= float(N**2)   # Mean magnetization
    # Compute derived values
    Cv = E_variance / T**2      # Speciffic Heat Capacitance
    chi = M_variance / T        # Magnetic Suceptibility

    return E_mean, M_mean, M_abs_mean, Cv, chi


def main():
    initspins = generateState(2)
    print(montecarlo(spins=initspins, T=1, trials=1e7))
    return


if __name__ == '__main__':
    main()
