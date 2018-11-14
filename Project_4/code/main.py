from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import metropolis as met
from scipy.constants import k


def figsetup(title, xlab, ylab, fname, legend=True, show=False, tightlayout=True):
    """
    Sets up and saves figure for usage in report
    usage:
    plot(...)
    plot(...)
    figsetup("filename")
    """
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(fname)
    plt.title(title)
    if tightlayout is True:
        plt.tight_layout()
    if legend is True:
        plt.legend()
    plt.savefig("../figs/" + fname + ".pdf")
    if show is False:
        plt.close()
    else:
        plt.show()
    return

def E_mean(Temp):
    beta = 1.0 / (k * Temp) 
    return -8 * np.sinh(8 * beta) / (np.cosh(8 * beta) + 3) 

def prob_1():
    

    return


def prob_2():

    L = 25
    n_temps = 10
    temps = np.linspace(2.2, 2.4, n_temps)
    initstate = met.generateState(L, 1)
    E = np.zeros(n_temps)
    E2 = np.zeros(n_temps)
    for i in range(n_temps):
        E[i], E2[i] = met.montecarlo(spins=initstate, T=temps[i], trials=1e6)

    plt.plot(temps, E / L**2)
    plt.show()

    return

def main():

    prob_2()
    return



if __name__ == '__main__':
    main()