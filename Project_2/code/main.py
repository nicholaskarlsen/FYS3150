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
import test_functions as tests


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


def sort_eigenpair(in_vals, in_vecs):
    "Sorts eigenpair such that eigenvals are in ascending order"
    out_vals = np.copy(in_vals)
    out_vecs = np.copy(in_vecs)

    permute = out_vals.argsort()
    out_vals = out_vals[permute]
    out_vecs = out_vecs[:, permute]

    return out_vals, out_vecs


def potential1(x):
    "Used for when there is no potential."
    return 0


def construct_matrix(dim, varMax=1.0, potential=potential1):
    step = varMax / (dim + 1)    # h, step size
    num = np.array(range(1, dim + 1))         # i = 1, ..., N
    var = np.zeros(dim) + num * step        # rho = rho_0 + ih

    d = 2.0 / step**2 * np.ones(dim) + potential(var)
    a = - 1.0 / step**2 * np.ones(dim)
    output = scipy.sparse.spdiags([a, d, a], [-1, 0, 1], dim, dim).toarray()  # Generates matrix

    return output, var


def potential2(x):
    "potential in question (d)"
    return x * x


def exercise_c():
    "performs jacobi for no potential, and times it"
    list_of_N = [5, 10, 20, 50, 70, 100, 150, 200, 300, 400, 500, 750, 1000]
    no_itterations = []
    time_taken = []
    for N in list_of_N:
        A, rho = construct_matrix(dim=N)
        A_eval, A_evec, counter, time = je.jacobi_solve(A, diagnostics=True)
        no_itterations.append(counter)
        time_taken.append(time)
    # Because arrays are easier to deal with
    no_itterations = np.array(no_itterations)
    time_taken = np.array(time_taken)
    list_of_N = np.array(list_of_N)

    plt.figure(figsize=[5, 5])
    plt.semilogy(list_of_N, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam", xlab="N",
             ylab="Time elapsed [s]", fname="q2c_time")

    plt.figure(figsize=[5, 5])
    plt.semilogy(list_of_N, no_itterations, "o--")
    figsetup(title="No. Itterations for solution of N-step Buckling beam", xlab="N",
             ylab="No. Itterations", fname="q2c_count")

    plt.figure(figsize=[5, 5])
    plt.loglog(list_of_N, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam", xlab="N",
             ylab="Time elapsed [s]", fname="q2c_timeloglog")

    plt.figure(figsize=[5, 5])
    plt.loglog(list_of_N, no_itterations, "o--")
    figsetup(title="No. Itterations for solution of N-step Buckling beam", xlab="N",
             ylab="No. Itterations", fname="q2c_countloglog")

    plt.figure(figsize=[5, 5])
    plt.plot(list_of_N, time_taken, "o--")
    figsetup(title="Time taken for solution of N-step Buckling beam", xlab="N",
             ylab="Time elapsed [s]", fname="q2c_timenormal")

    plt.figure(figsize=[5, 5])
    plt.plot(list_of_N, no_itterations, "o--")
    figsetup(title="No. Itterations for solution of N-step Buckling beam", xlab="N",
             ylab="No. Itterations", fname="q2c_countnormal")

    return


def exercise_d():
    "calls pertaining to question d"
    print "- Solving for question d"
    N = 1000
    A, rho = construct_matrix(dim=N, varMax=6, potential=potential2)
    A_eval, A_evec = je.jacobi_solve(A)

    A_eval, A_evec = sort_eigenpair(A_eval, A_evec)

    plt.figure(figsize=[5, 5])

    for n in [0, 1, 2]:
        plt.plot(rho, A_evec[:, n], label="$\\lambda=$%.4f" % A_eval[n])

    figsetup(title="Dimensionless wavefunction for first 3 eigenstates", xlab="$\\rho$", ylab="$u(\\rho)$",
             fname="question2d%i" % N)
    print A_eval[:3]
    print "- Done Solving for question d"
    return


def exercise_e():
    """
    calls pertaining to question e
    NOTE: when writing construct_matrix(), i didnt plan on the potentials needing
    multiple inputs. So when needing to solve the system for different omega_r,
    i instead opted to define a new potential3(x) function within the loop
    in which i am solving the system, so that i can update the _omega_r variable as a
    'global variable' rather than feeding it to the function. Whilst this would probably
    get me shot over at the informatics building, i'd rather that, than have to restructure
    this code another time.
    """
    print "- Solving for question e"
    N = 500
    plt.figure(figsize=[5, 5])
    for _omega_r in [0.01, 0.5, 1, 5]:
        def potential3(x):  # If this confuses you (which it should), read the multiline comment  above
            o = _omega_r
            return o * o * x * x + 1.0 / x
        A, rho = construct_matrix(dim=N, varMax=10, potential=potential3)
        A_eval, A_evec = je.jacobi_solve(A)
        A_eval, A_evec = sort_eigenpair(A_eval, A_evec)
        plt.plot(rho, A_evec[:, 0], label="$\\omega_r=$%.2f" % _omega_r)

    figsetup(title="Dimensionless wavefunction for first eigenstates", xlab="$\\rho$", ylab="$u(\\rho)$",
             fname="question2e")
    print "- Done Solving for question e"
    return


if __name__ == '__main__':
    tests.run_tests()
    exercise_c()
    exercise_d()
    exercise_e()
    print "Done."
