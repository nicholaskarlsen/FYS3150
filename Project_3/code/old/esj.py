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
from horizons import *

#m = np.array([3.00348959632E-6, 1, 954.79194E-6])

planets = {"Sun": 10, "Mercury": 199, "Venus": 299, "Earth": 399}

x0, v0, m = fetch_data(planets)
"""
x0 = np.array([
    [0.8691953608497864, 0.4844146835619565, -2.221389262494324e-05],
    [0, 0, 0],
    [-2.602835339639919, -4.695935682401297, 0.07774370082326859]
])

v0 = np.array([
    [-0.008654353274152921, 0.01495666860203449, -4.070416896990818e-07],
    [0, 0, 0],
    [0.006514074726715846, -0.003306388774089694, -0.0001320879214417232]
]) * 365.25
"""

N = 1e6
tn = 1
#names = ["Earth", "Sun", "Jupiter"]
names = ["Sun", "Mercury", "Venus", "Earth"]

esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
esys.velocityverlet(esys.gravity)
esys.plot("x")