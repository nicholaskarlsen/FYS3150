# Python v2.7
# Solves the SIRS model as a set of ODE's using RK4
# By Nicholas Karlsen

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from numba import njit


def RK4(y, h, func):
    k1 = h * func(t, y)
    k2 = h * func(t + h / 2.0, y + k1 / 2.0)
    k3 = h * func(t + h / 2, y + k2 / 2.0)
    k4 = h * func(t + h, y + k3)

    y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    return y_next


def ODE_SIRS(S0, I0, R0, a, b, c, N, tN):
    """
    S0 : Initial suceptibles
    I0 : Initial Infected
    R0 : Initial recovered
    N  : No. integration points
    tN : Simulation Time
    """
    init_population = S0 + I0 + R0

    t = np.linspace(0, tN, N)   # Time
    h = tN / N                  # Step Size
    N = int(N)                  # Ensure Int

    # Initialize Arrays
    S = np.zeros(N)  # Suceptible
    I = np.zeros(N)  # Infected
    R = np.zeros(N)  # Recovered
    population = np.zeros(N)
    S[0] = S0
    I[0] = I0
    R[0] = R0
    population[0] = sum([S[0], I[0], R[0]])

    def dSdt(i):
        return c * R[i] - a * S[i] * I[i] / population[i]

    def dIdt(i): return a * S[i] * I[i] / population[i] - b * I[i]

    def dRdt(i): return b * I[i] - c * R[i]

    for i in range(N - 1):
        population[i + 1] = sum([S[i], I[i], R[i]])
        S[i + 1] = S[i] + h * dSdt(i)
        I[i + 1] = I[i] + h * dIdt(i)
        R[i + 1] = R[i] + h * dRdt(i)

    return np.array([S, I, R, t])


def main():
    params = ["S", "I", "R"]
    plt.figure(1, figsize=[8, 8])
    for b in [1, 2, 3, 4]:
        SIRt = ODE_SIRS(S0=300, I0=100, R0=0, a=4, b=b, c=0.5, N=1e3, tN=10)
        for param in [0, 1, 2]:
            plt.subplot(2, 2, b)
            plt.plot(SIRt[-1], SIRt[param], label=params[param])
            plt.legend(loc="center right")
            plt.xlabel("Time")
            plt.ylabel("No. people")
            plt.title("b = %i" % b)
    plt.axis("tight")
    plt.show()

    return


if __name__ == '__main__':
    main()
