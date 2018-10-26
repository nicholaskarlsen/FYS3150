from __future__ import division
#from horizons import *
from get_initconds import *
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits import mplot3d
from n_solver import *
import time
import os.path


def figsetup(title, xlab, ylab, fname, legend=True, show=False):
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
    if legend is True:
        plt.legend()
    plt.savefig("../figs/" + fname + ".pdf")
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def ex_a():
    return


def ex_b():
    x0 = np.array([[1, 0]], dtype=np.float64)             # [AU]
    v0 = np.array([[0, 2 * np.pi]], dtype=np.float64)     # [AU/yr]
    m = np.array([3.00348959632E-6], dtype=np.float64)
    names = ['earth']
    list_of_N = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    timing = [[], [], []]

    N_exponent = 0
    method_name = ["eulercromer", "velocityverlet", "eulerforward"]  # for filename purposes.
    for N in list_of_N:
        N_exponent += 1
        for i in range(3):  # Methods numbered 0, 1, 2 corresponding to list above
            # Initialize the solver
            # Not strictly required to re-initialize each time as initial conds. arent altered
            inst = n_solver(initPos=x0, initVel=v0, mass=m, N=N, tn=10)

            t_start = time.time()                    # Start timing algorithm
            inst.solarsystem(method=i, system=1)     # Solve using solarsystem model for method=i
            timing[i].append(time.time() - t_start)  # Record time

            pos, vel = inst.get()       # Fetch arrays

            # Arrays to store speed & radius values for each t
            r = np.zeros(len(pos[0]))
            s = np.zeros(len(vel[0]))
            # Calculate r, s for each t
            for j in range(len(pos[0])):
                r[j] = np.linalg.norm(pos[0, j])
                s[j] = np.linalg.norm(vel[0, j])

            plt.figure(figsize=[2.5, 2.5])
            for j in xrange(len(x0)):  # Plot trajectories
                plt.plot(pos[j, :, 0], pos[j, :, 1], label=names[j])
            plt.axis("equal")
            figsetup(title="N = $10^%i$" % N_exponent, xlab="x [AU]", ylab="y [AU]",
                     fname="ex_b_orbit_%s_%i" % (method_name[i], N_exponent), legend=False)

            plt.figure(figsize=[5, 3])
            plt.plot(s)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            figsetup(title="Speed of earth", xlab="Integration point, N", ylab="|v| [Au/Yr]", fname="ex_b_speed_%s_%i" % (method_name[i], N_exponent))

            plt.figure(figsize=[5, 3])
            plt.plot(r)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            figsetup(title="Distance between earth & sun", xlab="Integration point, N",
                     ylab="|r| [Au]", fname="ex_b_radius_%s_%i" % (method_name[i], N_exponent))

            plt.figure(figsize=[5, 3])
            plt.plot(r * s * m[0])
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            figsetup(title="Angular momentum of earth", xlab="Integration point, N",
                     ylab="L [$M_\\bigodot$Au^2/Yr]", fname="ex_b_angularmomentum_%s_%i" % (method_name[i], N_exponent))
    plt.figure(figsize=[5, 5])
    for j in [1, 2]:
        plt.loglog(list_of_N, timing[j], "o--", label=method_name[j])
    figsetup(title="Timing Algorithms", xlab="N", ylab="Time [s]", fname="timing_earthsun")

    return


def ex_c1():
    x0 = np.array([[1, 0]], dtype=np.float64)             # [AU]
    m = np.array([3.00348959632E-6], dtype=np.float64)
    names = ['earth']
    for speed in [1, 2, 3, 3.5, 4]:
        v0 = np.array([[0, speed * np.pi]], dtype=np.float64)     # [AU/Yr]
        inst = n_solver(initPos=x0, initVel=v0, mass=m, N=1e5, tn=10)
        inst.solarsystem(method=1, system=1)  # Velocity verlet & fixed sun
        pos, vel = inst.get()
        plt.plot(pos[0, :, 0], pos[0, :, 1], label="$|v_0|=%.2f \\pi$" % speed)
    plt.close()
    return


def ex_c2():
    x0 = np.array([[1, 0]], dtype=np.float64)             # [AU]
    m = np.array([3.00348959632E-6], dtype=np.float64)
    names = ['earth']
    for b in [2, 2.25, 2.5, 2.75, 3]:
        v0 = np.array([[0, 2 * np.pi]], dtype=np.float64)     # [AU/Yr]
        inst = n_solver(initPos=x0, initVel=v0, mass=m, N=1e5, tn=10, beta=b)
        inst.solarsystem(method=1, system=1)  # Velocity verlet & fixed sun
        pos, vel = inst.get()
        plt.plot(pos[0, :, 0], pos[0, :, 1], label="$\\beta = %.1f$" % b)
    plt.show()
    return


def ex_d():
    m = np.array([1, 3.00348959632E-6], dtype=np.float64)

    lst = [10, 399]

    names = ['sun', 'earth']
    x0, v0 = get_data(lst)
    N = 1e5
    tn = 1

    inst = n_solver(initPos=x0, initVel=v0, mass=m, N=N, tn=tn)
    inst.solarsystem()
    pos, vel = inst.get()
    ax = plt.axes(projection='3d')
    for i in range(len(lst)):
        xline = pos[i, :, 0]
        yline = pos[i, :, 1]
        zline = pos[i, :, 2]
        ax.plot3D(xline, yline, zline, label=names[i])
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU]')
    plt.legend()
    plt.savefig("../figs/exd_var_beta.pdf")
    plt.show()
    return


def ex_e():
    lst = [399, 599]
    names = ["Earth", "jupiter"]
    x0, v0 = get_data(lst)
    for yrs in [5, 10]:#, 15, 20, 25, 30]:
        for multi in [1, 10, 1000]:
            m = np.array([3.00348959632E-6, multi * 954.79194E-6], dtype=np.float64)
            inst = n_solver(initPos=x0, initVel=v0, mass=m, N=1e6, tn=yrs)
            inst.solarsystem(method=1, system=1)  # Velocity verlet & fixed sun
            pos, vel = inst.get()
            ax = plt.axes(projection='3d')
            for i in range(len(lst)):
                xline = pos[i, :, 0]
                yline = pos[i, :, 1]
                zline = pos[i, :, 2]
                ax.plot3D(xline, yline, zline, label=names[i])
            ax.set_xlabel('x [AU]')
            ax.set_ylabel('y [AU]')
            ax.set_zlabel('z [AU]')
            plt.legend()
            plt.savefig("../figs/exe_earth_jupiter_%i_%iyrs.pdf" % (multi, yrs))
            plt.close()

    return


def ex_f():
    m = np.array([1, 0.16601E-6, 2.4478383E-6,
                  3.00348959632E-6, 0.3227151E-6, 954.79194E-6,
                  285.8860E-6, 43.66244E-6, 51.51389E-6,
                  0.007396E-6], dtype=np.float64)

    lst = [10, 199, 299, 399, 499,
           599, 699, 799, 899]

    names = ['sun', 'mercury', 'venus', 'earth', 'mars', 'jupiter',
             'satur', 'uranus', 'neptune']
    x0, v0 = get_data(lst, referenceFrame="500@0")
    N = 1e6
    tn = 165

    inst = n_solver(initPos=x0, initVel=v0, mass=m, N=N, tn=tn)
    inst.solarsystem()
    pos, vel = inst.get()
    ax = plt.axes(projection='3d')
    for i in range(len(lst)):
        xline = pos[i, :, 0]
        yline = pos[i, :, 1]
        zline = pos[i, :, 2]
        ax.plot3D(xline, yline, zline, label=names[i])
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU]')
    plt.legend()
    plt.savefig("../figs/all_planets3d.pdf")
    plt.show()
    return


def ex_g():
    N = int(1e8)   # takes 7131.14 s
    tn = 100
    if os.path.isfile('mercurypos.npy') is False:
        m = np.array([0.16601E-6], dtype=np.float64)
        names = ['mercury']
        x0 = [[0.307499, 0]]
        v0 = [[0, 12.44]]
        inst = n_solver(initPos=x0, initVel=v0, mass=m, N=N, tn=tn)
        t0 = time.time()
        inst.solarsystem(method=1, system=2)
        print "Time taken: %.2f s" % (time.time() - t0)
        pos, vel = inst.get()
        np.save('mercurypos.npy', pos)
        np.save('mercuryvel.npy', vel)
    if os.path.isfile('mercury_rad.npy') is False:
        pos = np.load('mercurypos.npy')
        r = np.zeros(N)
        for i in xrange(N):
            r[i] = np.linalg.norm(pos[0, i])
        np.save('mercury_rad.npy', r)

    r = np.load('mercury_rad.npy')
    i = int(N - 2)
    while r[i] < r[i + 1]:
        i -= 1
    print i

    plt.figure(figsize=[5, 5])
    plt.plot(r)
    plt.plot(i, r[i], 'rx')
    plt.xlim(i - 100000, N)
    figsetup(title="Distance between mercury and sun", xlab="Integration point, N", ylab="|$\\vec r$| [Au/Yr]", fname="ex_g_radius")

    pos = np.load('mercurypos.npy')
    print np.arctan(pos[0, i, 1] / pos[0, i, 0])

    return


def main():
    # ex_a()
    # ex_b()
    # ex_c1()
    # ex_c2()
    # ex_d()
    #ex_e()
    #ex_f()
    ex_g()

    return


if __name__ == '__main__':
    main()
