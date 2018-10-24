from __future__ import division, print_function
import numpy as np
from numba import jitclass          # import the decorator
from numba import int32, float32, float64, uint8    # import the types

# Declaring types for @jitclass
spec = [
    ('initPos', float64[:, :]),
    ('initVel', float64[:, :]),
    ('mass', float64[:]),
    ('N', int32),
    ('tn', float32),
    ('h', float64),
    ('numBodies', int32),
    ('dim', int32),
    ('G', float64),
    ('pos', float64[:, :, :]),
    ('vel', float64[:, :, :]),
    ('method', uint8)
]


#@jitclass(spec)
class n_solver(object):
    def __init__(self, initPos, initVel, mass, N, tn, method=1):
        # Input variables
        self.initPos = initPos      # Initial condition [AU]
        self.initVel = initVel      # Initial velocity [Au/Yr]
        self.mass = mass            # Array of masses [Solar mass]
        self.N = int(N)                  # Number of integration points
        self.tn = tn                # Simulation time [Yr]
        self.h = self.tn / self.N
        # Useful constants
        self.numBodies = len(self.initPos)
        self.dim = len(self.initPos[0])  # Number of dimensions, should be 2 or 3
        self.G = 4 * np.pi**2           # Gravitational Constant [AU^3 yr^-2 M_sun^-1]

        # Using .self variables in the shape arg causes error with jitclass
        self.pos = np.zeros(shape=(len(initPos), int(N), len(initPos[0])), dtype=np.float64)
        self.vel = np.zeros(shape=(len(initPos), int(N), len(initPos[0])), dtype=np.float64)
        self.pos[:, 0] = self.initPos
        self.vel[:, 0] = self.initVel

    def get(self):
        print(self.method)
        return self.pos, self.vel, self.mass


if __name__ == '__main__':
    m = np.array([1, 3.00348959632E-6], dtype=np.float64)
    x0 = np.array([[0, 0], [1, 0]], dtype=np.float64)
    v0 = np.array([[0, 0], [0, 2 * np.pi]], dtype=np.float64)
    names = ["Sun", "Earth"]
    #planets = {"Earth": 399, "Sun": 10}
    #x01, v01, m1 = hori.fetch_data(jpl_id=planets)
    N = 1e5
    tn = 1

    esys = n_solver(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, method="ec")
    pos, vel, mass = esys.get()
