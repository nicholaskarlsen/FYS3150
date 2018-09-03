import numpy as np
import matplotlib.pyplot as plt

# Initial Conditions & parameters
n = 1e3
h = 1.0 / float(n)
u0 = 0  # u(0) = 0
u1 = 0  # u(1) = 0
x = np.linspace(0, 1, h)


def analyticSolution(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10*x)

def tridiagonal(n, dNum=2, eNum=-1):
    """
    Generates the diagonal elements of a tridiagonal [nxn] matrix
    with identical elements along the diagonal as two vectors.
    """
    if n < 4:  # Matrix does not make sense for n>4
        raise ValueError("n too small")

    d = np.ones(n) * dNum
    e = np.ones(n - 2) * eNum

    return [d, e]

xvals = np.linspace(0, 1, 1e3)

plt.plot(xvals, analyticSolution(xvals), label="Analytic Solution")
plt.legend()
plt.show()