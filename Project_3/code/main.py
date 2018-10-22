#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15
from __future__ import division  # Nobody expects the integer division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from solarsystem import *

# Contains calls which generate all figs & data for report.


def figsetup(title, xlab, ylab, fname, show=False):
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
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
    for N in tn * np.array([1e1, 1e2, 1e3, 1e4, 1e5]):
        esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
        esys.velocityverlet(esys.gravity)
        esys.changeref(1)
        pos, vel = esys.get()
        epos = pos[0]
        plt.plot(epos[:, 0], epos[:, 1], label="%.1e" % N)
    plt.legend()
    plt.xlim(-1.1, 1.1)
    plt.ylim(-1.1, 1.1)
    plt.show()   
    return


if __name__ == '__main__':
    ex_c()
