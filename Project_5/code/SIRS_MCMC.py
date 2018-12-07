from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from numba import njit


def genState(NumS, NumI, NumR):
    # generates shuffled array with population
    # entries; 'S', 'I', 'R'

    # Start by filling up an array with "S"
    initState = np.empty(NumS, dtype=str)
    initState.fill("S")
    # Repeat for a new array with "I"
    tmp = np.empty(NumI, dtype=str)
    tmp.fill("I")
    # Append to Previous array
    initState = np.append(initState, tmp)
    # Repeat for "R"
    tmp = np.empty(NumR, dtype=str)
    tmp.fill("R")
    initState = np.append(initState, tmp)
    # Shuffle the array randomly
    initState = np.random.permutation(initState)
    return initState


def countState(state):
    NumS = np.in1d(state, 'S').sum()
    NumI = np.in1d(state, 'I').sum()
    NumR = np.in1d(state, 'R').sum()
    return NumS, NumI, NumR

def main(S0, I0, R0, a, b, c, trials=2, stop_time=10):
    N = S0 + I0 + R0        # Initial population

    dt = np.min([           # Step size
        4.0 / (a * N),
        1.0 / (b * N),
        1.0 / (c * N)])

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


def testFunc(NumTests=1e3):
    """ 
    Test functionality of genState() & countState() by asserting that they
    should yield the same results
    """
    NumTests = int(NumTests)
    msg = "genStates & countState yield different results"
    for i in range(NumTests):
        x1 = np.random.randint(0, 1000)
        y1 = np.random.randint(0, 1000)
        z1 = np.random.randint(0, 1000)
        tmpState = genState(x1, y1, z1)
        x2, y2, z2 = countState(tmpState)

        if x1 != x2:
            raise ValueError(msg)
        elif y1 != y2:
            raise ValueError(msg)
        elif z1 != z2:
            ValueError(msg)

    return


if __name__ == '__main__':
    t, S, I, R = main(S0=300, I0=100, R0=0, a=4, b=1, c=0.5, stop_time=10, trials=100)
    plt.plot(t, S, color="blue", alpha=.1)
    plt.plot(t, I, color="red", alpha=.1)
    plt.plot(t, R, color="green", alpha=.1)

    plt.legend()
    plt.show()

    testFunc(1)
