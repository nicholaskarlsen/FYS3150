#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15
from __future__ import division  # Nobody expects the integer division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

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

if __name__ == '__main__':
        
    fig = plt.figure()
    ax = plt.axes(projection='3d')#

    # Data for a three-dimensional line
    zline = np.linspace(0, 15, 1000)
    xline = np.sin(zline)
    yline = np.cos(zline)
    ax.plot3D(xline, yline, zline)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()