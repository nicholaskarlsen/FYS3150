import numpy as np


def generateState(N):
    # function wrapper which returns [N, N] array of random values in [-1, 1]
    array = np.random.choice([-1, 1], size=(N, N))
    return array
