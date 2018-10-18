# This script solves the 2-body problem using a simple approach in two dimensions
from __future__ import division
import numpy as np
from scipy import constants as const
from horizons import fetch_data
from solver import diffeqsolver
import matplotlib.pyplot as plt

G = const.G

x0 = np.array([const.au, 0])  # (x, y) [AU]
v0 = np.array([0, np.sqrt(4 * np.pi * const.au**2 / (362.25 * 24 * 3600))])


tn = 3600 * 24 * 365
N = int(1e6)
dt = tn / N

pos = np.zeros([N, 2])
vel = np.zeros([N, 2])

pos[0] = x0
vel[0] = v0


def gravity(r, m=6e24):
    return -G * m / r**2


def velocityverlet(current_pos, current_vel):
    a1 = gravity(current_pos)   # a_i
    next_pos = current_pos + current_vel * dt + 0.5 * a1 * dt**2
    a2 = gravity(next_pos)
    next_vel = current_vel + 0.5 * dt * (a1 + a2)
    return next_pos, next_vel


def eulercromer(current_pos, current_vel):
    next_vel = current_vel + dt * gravity(current_pos)
    next_pos = current_pos + dt * next_vel
    return next_pos, next_vel


for i in xrange(N - 1):
    for dim in xrange(2):
        pos[i + 1][dim], vel[i + 1][dim] = eulercromer(pos[i][dim], vel[i][dim])


plt.plot(pos[:][0], pos[:][1])
plt.show()
