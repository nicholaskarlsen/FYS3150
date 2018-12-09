# Python v2.7
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def main(S0, I0, R0, a, b, c, trials=2, stop_time=10):
    N = S0 + I0 + R0        # Initial population

    dt = np.min([4.0 / (a * N), 1.0 / (b * N), 1.0 / (c * N)])  # Step size

    t = np.arange(0, stop_time, step=dt)    # Time

    S = np.zeros([len(t), trials], dtype=int)  # Suceptibles
    I = np.zeros([len(t), trials], dtype=int)  # Infected
    R = np.zeros([len(t), trials], dtype=int)  # Recovered

    # Set initial conditions
    S[0, :] = S0
    I[0, :] = I0
    R[0, :] = R0

    for i in range(len(t) - 1):
        for j in range(trials):
            # Initialize next step
            S[i + 1, j] = S[i, j]
            I[i + 1, j] = I[i, j]
            R[i + 1, j] = R[i, j]

            # Transition probabilities
            S_I = (a * S[i, j] * I[i, j] * dt) / (S[i, j] + I[i, j] + R[i, j])   # P(S->I)
            I_R = b * I[i, j] * dt                                   # P(I->R)
            R_S = c * R[i, j] * dt                                   # P(R->S)

            # if lt <= S[i]:
            # (S -> I)
            if np.random.random() < S_I:
                S[i + 1, j] -= 1
                I[i + 1, j] += 1

            # elif lt > S[i] + I[i]:
                # (R -> S)
            if np.random.random() < R_S:
                R[i + 1, j] -= 1
                S[i + 1, j] += 1

            # else:
                # (I -> R)
            if np.random.random() < I_R:
                I[i + 1, j] -= 1
                R[i + 1, j] += 1

    return t, S, I, R


if __name__ == '__main__':
    t, S, I, R = main(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10, trials=100)
    plt.plot(t, S, color="blue", alpha=.1)
    plt.plot(t, I, color="red", alpha=.1)
    plt.plot(t, R, color="green", alpha=.1)

    plt.legend()
    plt.show()
