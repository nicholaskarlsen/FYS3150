import numpy as np
import matplotlib.pyplot as plt
from thomas_algorithm import *
from gausselim import *

# Initial Conditions & parameters
n = 1e3
h = 1.0 / float(n)
u0 = 0  # u(0) = 0
u1 = 0  # u(1) = 0
x = np.linspace(0, 1, h)

def f_func(x):
    """ Source term """
    return 100 * np.exp(-10 * x)


def analyticSolution(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


xvals = np.linspace(0, 1, 1e3)

plt.plot(xvals, analyticSolution(xvals), label="Analytic Solution")
for num in [10, 100, 1000]:
    plt.plot(gauss_general(num, f_func)[0], gauss_general(num, f_func)[1], "x-", label="Tridiagonal Solution (n=%i)" % num)
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.show()
