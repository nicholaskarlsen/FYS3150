# Avoid integer division + jit works better with print func.
from __future__ import division, print_function
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import rcParams
from horizons import *
rcParams.update({'figure.autolayout': True})

planets = {"Sun": 10, "Earth": 399}# "Venus": 299}
x0, v0, m = fetch_data(planets)

@jit(nopython=True)
def earthsun(initPos, initVel, mass, N, tn):
    h = tn / N
    N = int(N)                  # Number of integration points
    # Useful constants
    numBodies = len(initPos)
    dim = len(initPos[0])  # Number of dimensions, should be 2 or 3
    G = 4 * np.pi**2           # Gravitational Constant [AU^3 yr^-2 M_sun^-1]

    # Using .self variables in the shape arg causes error with jitclass
    pos = np.zeros(shape=(len(initPos), int(N), len(initPos[0])), dtype=np.float64)
    vel = np.zeros(shape=(len(initPos), int(N), len(initPos[0])), dtype=np.float64)

    # Change to sun frame of reference
    for i in range(len(pos)):
        initPos[i] -= initPos[0]
        initVel[i] -= initVel[0]

    pos[:, 0] = initPos
    vel[:, 0] = initVel

    for t in xrange(N - 1):
        acc =  G * pos[1][i] / np.linalg.norm(pos[i]) **3
        vel[1][i + 1] = vel[1][i] + h
        pos[1][i + 1] = pos[1][i] + h * vel[1][i + 1]


    return pos, vel, mass



pos, vel, mass = earthsun(x0, v0, m, 1e5, 1) 

ax = plt.axes(projection='3d')

xline = pos[1, :, 0]
yline = pos[1, :, 1]
zline = pos[1, :, 2]
ax.plot3D(xline, yline, zline)
ax.set_xlabel('x [AU]')
ax.set_ylabel('y [AU]')
ax.set_zlabel('z [AU]')
plt.legend()
plt.show()