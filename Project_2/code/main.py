# Generates all of the figures and data for the report
# Python Version: Python 2.7.15

"""
Tried to avoid reusing variable names as much as possible due to a catasrophic
bug in a previous version of these scripts. Therefore, there may not always be
a 1:1 correspondance between the names in the project text and the
variable names in this script.
"""

from __future__ import division  # Nobody expects the integer division

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.sparse
import jacobi_eigensolver as je


def figsetup(title, xlab, ylab, fname, show=False):
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
    plt.tight_layout()
    plt.title(title)
    plt.legend()
    plt.savefig("../figs/" + fname + ".png", dpi=250)
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def potential1(x):
    "Used for when there is no potential."
    return 0


def potential2(x):
    return x * x


def sort_eigenpair(in_vals, in_vecs):
    "Sorts eigenpair such that eigenvals are in ascending order"
    out_vals = np.copy(in_vals)
    out_vecs = np.copy(in_vecs)

    permute = out_vals.argsort()
    out_vals = out_vals[permute]
    out_vecs = out_vecs[:, permute]

    return out_vals, out_vecs


def construct_matrix(dim, varMax=1.0, potential=potential1):
    step = varMax / (dim + 1)    # h, step size
    num = np.array(range(1, dim + 1))         # i = 1, ..., N
    var = np.zeros(dim) + num * step        # rho = rho_0 + ih

    d = 2.0 / step**2 * np.ones(dim) + potential(var)
    a = - 1.0 / step**2 * np.ones(dim)
    output = scipy.sparse.spdiags([a, d, a], [-1, 0, 1], dim, dim).toarray()  # Generates matrix

    return output, var

def exercise_d():
    "calls pertaining to question d"
    print "--- Solving for question d ---"
    N = 1000
    A, rho = construct_matrix(dim=N, varMax=10, potential=potential2)
    A_eval, A_evec = je.jacobi_solve(A)

    A_eval, A_evec = sort_eigenpair(A_eval, A_evec)

    plt.figure(figsize=[5,5])
    
    for n in [0, 1, 2]:
        plt.plot(rho, A_evec[:, n], label="$\\lambda=$%.4f" % A_eval[n])
    
    figsetup(title="Dimensionless wavefunction for first 3 eigenstates", xlab="$\\rho$", ylab="$u(\\rho)$",
             fname="question2d")
    print A_eval[:3]
    print "--- Done Solving for question d ---"
    return


if __name__ == '__main__':
    #exercise_d()

    print "Done."