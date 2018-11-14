# Python Version 2.7.15
from __future__ import division, print_function  # Nobody expects the integer division
import numpy as np
from numba import jit


def generateState(N, spin=False):
    # function wrapper which returns [N, N] array of random values in [-1, 1]
    # unless the spin parameter is set, in which case it will return a filled array
    if spin is False:
        array = np.random.choice([-1, 1], size=(N, N))
    else:
        array = np.full(shape=(N, N), fill_value=spin)

    return array



@jit(nopython=True)
def montecarlo(spins, T, trials):
    """
    N: size of square lattice
    initState: Initial state of system, [N, N] array with entries +/-1
    T: Temperature of the system
    trials: Number of montecarlo cycles
    """

    def periodic(i, limit, add):
        """
        Choose correct matrix index with periodic
        boundary conditions

        Input:
        - i:     Base index
        - limit: Highest \"legal\" index
        - add:   Number to add or subtract from i
        """
        return (i + limit + add) % limit

    N = len(spins)                             # Ensure that N is int
    trials = int(trials)                    # ...           is int

    # If no init state is given
    #spins = np.ones((N, N), dtype=np.int8)  # Create initial state of spin-ups

    E = 0                                   # Compute its energy
    for j in xrange(N):
        for i in xrange(N):
            E -= spins[i, j] * (spins[periodic(i, N, -1), j] + spins[i, periodic(j, N, 1)])

    w = np.zeros(17, dtype=np.float64)
    for de in xrange(-8, 9, 4):  # include +8
        w[de + 8] = np.exp(-de / T)

    E2_mean = 0   # Expectation value of E^2
    E_mean = 0    # Expectation value of E

    # Monte carlo sim
    for i in xrange(trials):
        for s in xrange(int(N**2)):
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

        E_mean += E
        E2_mean += E**2

    E_mean /= float(trials)
    E2_mean /= float(trials)
    # Calculate variance and normalize to per-point and temp
    E_variance = (E2_mean - E_mean * E_mean) / float(N**2 * T**2)
    # Normalize returned averages to per-point
    E_mean /= float(N**2)

    return E_mean, E2_mean


def main():
    initspins = generateState(2)
    print(montecarlo(spins=initspins, T=1, trials=1e4))
    return


if __name__ == '__main__':
    main()
