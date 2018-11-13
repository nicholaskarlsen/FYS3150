# Python Version 2.7.15
from __future__ import division  # Nobody expects the integer division
import numpy as np


def generateState(N, spin=False):
    # function wrapper which returns [N, N] array of random values in [-1, 1]
    # unless the spin parameter is set, in which case it will return a filled array
    if spin is False:
        array = np.random.choice([-1, 1], size=(N, N))
    else:
        array = np.full(shape=(N, N), fill_value=spin)

    return array

def montecarlo(initState, T, trials):
    """
    initState: Initial state of system, [N, N] array with entries +/-1
    T: Temperature of the system
    trials: Number of montecarlo cycles
    """
    spins = initState  

    return


def main():

    return


if __name__ == '__main__':
    main()
