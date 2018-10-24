# Avoid integer division + jit works better with print func.
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from numba import jitclass, jit          # import the decorator
from numba import int32, float32, float64, uint8    # import the types
import time
from get_initconds import *

# Wanted to implement jitclass to speed up code, but am not able to
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
    ('c', float64),
    ('pos', float64[:, :, :]),
    ('vel', float64[:, :, :]),
    ('beta', float64),
]


#@jitclass(spec)
class n_solver(object):
    def __init__(self, initPos, initVel, mass, N, tn, beta=2.0):
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
        self.G = 4 * np.pi**2            # Gravitational Constant [AU^3 yr^-2 M_sun^-1]
        self.c = 63197.8
        self.beta = beta                 # Exponent in gravity method (1/r^beta)

        # Using .self variables in the shape arg causes error with jitclass
        self.pos = np.zeros(shape=(len(initPos), int(N), len(initPos[0])), dtype=np.float64)
        self.vel = np.zeros(shape=(len(initPos), int(N), len(initPos[0])), dtype=np.float64)
        for i in xrange(self.numBodies):
            self.pos[i, 0] = self.initPos[i]
            self.vel[i, 0] = self.initVel[i]

    def gravity_fixedsun(self, planetIndex, timeIndex):
        """Fixes the sun at the origin & calculates gravitational attraction to it, and
        any other planets which may be present in the system """
        accel = np.zeros(self.dim)
        accel -= self.pos[planetIndex, timeIndex] * self.G / np.linalg.norm(self.pos[planetIndex, timeIndex])**(1 + self.beta)
        for j in xrange(self.numBodies):
            if j != planetIndex:  # Planet doesnt act on itself
                relPos = self.pos[planetIndex, timeIndex] - self.pos[j, timeIndex]
                accel -= (relPos * self.G * self.mass[j]) / np.linalg.norm(relPos) ** (1 + self.beta)
            else:
                pass
        return accel

    def gravity(self, planetIndex, timeIndex):
        accel = np.zeros(self.dim)
        for j in xrange(self.numBodies):
            if j != planetIndex:  # Planet doesnt act on itself#
                relPos = self.pos[planetIndex, timeIndex] - self.pos[j, timeIndex]
                accel -= (relPos * self.G * self.mass[j]) / np.linalg.norm(relPos) ** (1 + self.beta)
            else:
                pass
        return accel

    def gravity_relativistic(self, planetIndex, timeIndex):
        """Fixes the sun at the origin & calculates gravitational attraction to it, and
        any other planets which may be present in the system + adds relativistic
        correction """
        accel = np.zeros(self.dim)  # Used to store values

        # Calculate gravity between sun & planet
        l_vec = np.cross(self.pos[planetIndex, timeIndex], self.vel[planetIndex, timeIndex])  # Angular momentum
        l = np.linalg.norm(l_vec)  # Magnitude of angular momentum
        relCorr = 1 + ((3 * l**2) / (self.pos[planetIndex, timeIndex]**2 * self.c**2))  # correction factor
        # Computes gravitational acceleration due to fixed sun
        accel -= relCorr * self.pos[planetIndex, timeIndex] * self.G / np.linalg.norm(self.pos[planetIndex, timeIndex])**(1 + self.beta)

        for j in xrange(self.numBodies):
            if j != planetIndex:  # Planet doesnt act on itself
                relPos = self.pos[planetIndex, timeIndex] - self.pos[j, timeIndex]
                relVel = self.vel[planetIndex, timeIndex] - self.vel[j, timeIndex]

                # Calculate the relativistic correction
                l_vec = np.cross(relPos, relVel)  # Angular momentum
                l = np.linalg.norm(l_vec)  # Magnitude of angular momentum
                relCorr = 1 + ((3 * l**2) / (self.pos[planetIndex, timeIndex]**2 * self.c**2))  # correction factor

                accel -= relCorr * (relPos * self.G * self.mass[j]) / np.linalg.norm(relPos) ** (1 + self.beta)
            else:
                pass
        return accel

    def eulerforward(self, diffeq):
        """ 
        Solves N coupled differential equations using Forward Euler
        for a given differential equationi, diffeq
        diffeq: function(t, n)
        """
        print("Solving with Forward Euler algorithm with N=", self.N)
        for t in xrange(self.N - 1):
            for n in xrange(self.numBodies):
                acc = diffeq(n, t)
                self.vel[n][t + 1] = self.vel[n][t] + self.h * acc
                self.pos[n][t + 1] = self.pos[n][t] + self.h * self.vel[n][t]
        return

    def eulercromer(self, diffeq):
        """ 
        Solves N coupled differential equations using Euler-Cromer
        for a given differential equation
        """
        print("Solving with Euler-Cromer algorithm with N=", self.N)
        for t in xrange(self.N - 1):
            for n in xrange(self.numBodies):
                acc = diffeq(n, t)
                self.vel[n][t + 1] = self.vel[n][t] + self.h * acc
                self.pos[n][t + 1] = self.pos[n][t] + self.h * self.vel[n][t + 1]
        return

    def velocityverlet(self, diffeq):
        print("Solving with Velocity-Verlet algorithm with N=", self.N)
        for t in xrange(self.N - 1):
            a1 = np.zeros([self.numBodies, self.dim])
            a2 = np.zeros([self.numBodies, self.dim])
            for n in xrange(self.numBodies):  # First find pos+1 for all planets
                a1[n] = diffeq(n, t)
                self.pos[n][t + 1] = self.pos[n][t] + self.vel[n][t] * self.h + 0.5 * a1[n] * self.h ** 2
            for n in xrange(self.numBodies):  # Then find vel+1 for all planets
                a2[n] = diffeq(n, t + 1)
                self.vel[n][t + 1] = self.vel[n][t] + 0.5 * self.h * (a2[n] + a1[n])
        return

    def solarsystem(self, method=1, system=0):
        """
        Because of jitclass limiting nature ive had to resort to using a function
        wrapper that calls the code in a far less general manner than i originally intended.
        Further, JITclass doesnt allow input/output of strings (for whatever reason)
        hence "solvers" list.
        NOTE: JIT did not like this way of doing things either. wrapper remains
              because it makes calling the code easy.
        """
        diffeqs = [self.gravity, self.gravity_fixedsun, self.gravity_relativistic]
        solvers = [self.eulercromer, self.velocityverlet, self.eulerforward]
        solvers[method](diffeqs[system])
        return

    def get(self):
        # Used to fetch the data
        return self.pos, self.vel

    def plot(self):
        # Basic plotting method that works without input for both 2&3 dimensions
        if self.dim == 3:
            ax = plt.axes(projection='3d')
            for i in range(self.numBodies):
                # Data for a three-dimensional line
                xline = self.pos[i, :, 0]
                yline = self.pos[i, :, 1]
                zline = self.pos[i, :, 2]
                ax.plot3D(xline, yline, zline)
            ax.set_xlabel('x [AU]')
            ax.set_ylabel('y [AU]')
            ax.set_zlabel('z [AU]')
            plt.show()

        if self.dim == 2:
            for i in range(self.numBodies):
                plt.plot(self.pos[i, :, 0], self.pos[i, :, 1])
            plt.xlabel("x [AU]")
            plt.ylabel("y [AU]")
            plt.show()

        return
