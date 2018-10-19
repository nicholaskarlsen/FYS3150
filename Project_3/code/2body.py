# This script solves the 2-body problem using a simple approach in two dimensions
from __future__ import division
import numpy as np
from scipy import constants as const
from horizons import fetch_data
from solver import diffeqsolver
import matplotlib.pyplot as plt
from numba import jit

G = 4 * np.pi ** 2    # [AU^3 / (Yr^2 * M_sun)]

x0 = np.array([1, 0])           # (x, y) [AU]
v0 = np.array([0, 2 * np.pi])   # (vx, vy) [AU / yr]


tn = 1          # Timespan
N = int(1e6)    # Number of integration points
h = tn / N      # Step size

pos = np.zeros([N, 2])
vel = np.zeros([N, 2])

pos[0] = x0
vel[0] = v0


def gravity(pos, mass):
    return pos * G * mass / np.linalg.norm(pos)**3


def velocityverlet(diffeq):
    a1 = diffeq(pos[i])   # a_i
    pos[i + 1] = pos[i] + vel[i] * h + 0.5 * a1 * h**2
    a2 = diffeq(pos[i + 1])
    vel[i + 1] = vel[i] + 0.5 * h * (a1 + a2)
    return next_pos, next_vel


def eulercromer(diffeq):
    for i in xrange(N - 1):
        acc = diffeq(pos[i], 1)
        vel[i + 1] = vel[i] + h * acc
        pos[i + 1] = pos[i] + h * vel[i + 1]
    return pos, vel


for i in range(N - 1):
    r = np.linalg.norm(pos[i])
    for dim in [0, 1]:
        a = -pos[i][dim] * G / r ** 3
        vel[i + 1][dim] = vel[i][dim] + h * a
        pos[i + 1][dim] = pos[i][dim] + h * vel[i + 1][dim]

plt.plot(pos[:][0])
plt.plot(pos[:][1])
plt.show()
plt.plot(pos[:][0], pos[:][1])
plt.show()
