from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from numba import njit


def main(S0, I0, R0, a, b, c, stop_time=10):
    N = S0 + I0 + R0                                            # Initial population
    dt = np.min([4.0 / (a * N), 1.0 / (b * N), 1.0 / (c * N)])  # Step size
    t = np.arange(0, stop_time, step=dt)                        # Time

    S = np.zeros(len(t), dtype=int)  # Suceptibles
    I = np.zeros(len(t), dtype=int)  # Infected
    R = np.zeros(len(t), dtype=int)  # Recovered

    # Set initial conditions
    S[0] = S0
    I[0] = I0
    R[0] = R0

    for i in range(len(t) - 1):
        # Initialize next step
        S[i + 1] = S[i]
        I[i + 1] = I[i]
        R[i + 1] = R[i]

        # Probabilities
        S_I = a * S[i] * I[i] * dt  # P(S->I)
        I_R = b * I[i] * dt         # P(I->R)
        R_S = c * R[i] * dt         # P(R->S)

        lt = np.random.randint(S[i] + I[i] + R[i])  # Choses S, I or R
        r = np.random.random()                     # acceptance

        if lt <= S[i]:
            # (S -> I)
            if r < S_I:
                S[i + 1] -= 1
                I[i + 1] += 1

        elif lt >= S[i] + I[i]:
            # (R -> S)
            if r < R_S:
                R[i + 1] -= 1
                S[i + 1] += 1

        else:
            # (I -> R)
            if r < I_R:
                I[i + 1] -= 1
                R[i + 1] += 1

    return t, S, I, R


if __name__ == '__main__':
    t, S, I, R = main(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10)
    plt.plot(t, S, label="S")
    plt.plot(t, I, label="I")
    plt.plot(t, R, label="R")

    plt.legend()
    plt.show()
