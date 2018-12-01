# Python v2.7
# Solves the SIRS model as a set of ODEs using Runge-Kutta 4
# By Nicholas Karlsen

from __future__ import division, print_function  # Avoid unexpected integer divisions & printf for JIT
import numpy as np
from numba import njit


def main(S0, I0, R0, a, b, c, N, tN):

    t = np.linspace(0, tN, N)

    S = np.zeros(N)  # Suceptible
    I = np.zeros(N)  # Infected
    R = np.zeros(N)  # Recovered
    population = np.zeros(N)
    S[0] = S0
    I[0] = I0
    R[0] = R0
    population[0] = sum([S[0], I[0], R[0]])

    # Define

    def diffEqs(t, S, I, R):
        """
        Compute differentials for the SIRS system at a time t. Kept as a nested function
        to access the local namespace for a, b, c constants & allow for the usage of
        numba.jit.

        Parameters
        ----------
        t: Time
        S: Suceptibles
        I: Infected
        R: Recovered

        Returns
        -------
        Numpy array
            Array filled with differentials governing the system at time t.
        """
        dSdt = c * R - a * S * I / (S + I + R)
        dIdt = a * S * I / (S + I + R) - b * I
        dRdt = b * I - c * R

        return np.array([dSdt, dIdt, dRdt])

         

        return y_next
