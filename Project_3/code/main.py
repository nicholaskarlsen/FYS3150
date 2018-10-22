#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15
from __future__ import division  # Nobody expects the integer division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from solarsystem import *
from matplotlib import rcParams
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf

conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'
rcParams.update({'figure.autolayout': True})

# Contains calls which generate all figs & data for report.


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
    plt.savefig("../figs/" + fname + ".pdf")
    if show is False:
        plt.close()
    else:
        plt.show()
    return


def ex_c():
    "Calls generating figs for ex c"
    m = np.array([3.00348959632E-6, 1])
    x0 = np.array([[1, 0], [0, 0]])
    v0 = np.array([[0, 2 * np.pi], [0, 0]])
    #planets = {"Earth": 399, "Sun": 10}
    #x01, v01, m1 = hori.fetch_data(jpl_id=planets)
    names = ["Earth", "Sun"]
    tn = 10
    list_of_N = tn * np.array([1e1, 1e2, 1e3, 1e4, 1e5])
    for N in list_of_N:
        esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
        esys.velocityverlet(esys.gravity)
        esys.changeref(1)
        pos, vel = esys.get()
        epos = pos[0]
        plt.plot(epos[:, 0], epos[:, 1], label="%.1e" % N)

    def excplots(fname):
        plt.legend()
        plt.xlim(-1.1, 1.1)
        plt.ylim(-1.1, 1.1)
        plt.savefig("../figs/ex_c_" + fname + "_orbit.pdf")
        plt.close()
        # System is still "solved" for the last N
        pos, vel = esys.get()  # Returns the pos and vel arrays from the object
        epos, evel = pos[0], vel[0]
        Ek = np.zeros(int(list_of_N[-1]))   # kinetic energy
        Ep = np.zeros(int(list_of_N[-1]))   # potential energy
        _G = 4 * np.pi ** 2
        for i in xrange(int(list_of_N[-1])):
            r = np.linalg.norm(epos[i])
            v = np.linalg.norm(evel[i])
            Ek[i] = 0.5 * m[0] + v**2
            Ep[i] = - m[0] * _G / r

        plt.plot(Ek + Ep)
        plt.xlabel("i")
        plt.ylabel("Energy")
        plt.legend()
        plt.savefig("../figs/ex_c_" + fname + "_energy.pdf")
        plt.close()

    excplots("vv")
    for N in list_of_N:
        esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
        esys.eulercromer(esys.gravity)
        esys.changeref(1)
        pos, vel = esys.get()
        epos = pos[0]
        plt.plot(epos[:, 0], epos[:, 1], label="%.1e" % N)
    excplots("ec")

    return


def ex_d():
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    m = np.array([3.00348959632E-6, 1])
    x0 = np.array([[1, 0, 0], [0, 0, 0]])
    for speed in [0.5 * np.pi, np.pi, 2 * np.pi, 3 * np.pi]:
        v0 = np.array([[0, 2 * np.pi, speed], [0, 0, 0]])
        #planets = {"Earth": 399, "Sun": 10}
        #x01, v01, m1 = hori.fetch_data(jpl_id=planets)
        names = ["Earth", "Sun"]
        tn = 5
        N = tn * 1e5
        esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
        esys.velocityverlet(esys.gravity)
        pos, vel = esys.get()
        # Data for a three-dimensional line
        xline = pos[0, :, 0]
        yline = pos[0, :, 1]
        zline = pos[0, :, 2]
        ax.plot3D(xline, yline, zline, label="$v_0=%.3f$" % np.sqrt(speed**2 + 4 * np.pi**2))
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU]')
    plt.legend()
    plt.savefig("../figs/" + "escapevel" + ".pdf")
    plt.show()

    return


def ex_e():
    m = np.array([3.00348959632E-6, 1, 954.79194E-6])

    sun = Horizons(id="10", id_type="id", location="500@10")
    earth = Horizons(id="399", id_type="id", location="500@10")
    jupiter = Horizons(id="599", id_type="id", location="500@10")

    x0 = np.array([
        [0.8691953608497864, 0.4844146835619565, -2.221389262494324e-05],
        [0, 0, 0],
        [-2.602835339639919, -4.695935682401297, 0.07774370082326859]
    ])

    v0 = np.array([
        [-0.008654353274152921, 0.01495666860203449, -4.070416896990818e-07],
        [0, 0, 0],
        [0.006514074726715846, -0.003306388774089694, -0.0001320879214417232]
    ]) / (365.25)

    N = 1e6
    tn = 10
    names = ["Earth", "Sun", "Jupiter"]

    esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
    esys.velocityverlet(esys.gravity)
    esys.plot("earthsun")

    return


if __name__ == '__main__':
    # ex_c()
    # ex_d()
    ex_e()
