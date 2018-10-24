from __future__ import division  # Nobody expects the integer division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import horizons as hori

from numba import jitclass          # import the decorator
from numba import int32, float32, float64, char    # import the types

spec = [
    ('mass', float64[:]),
    ("N", int32),
    ("tn", float32),
    ("names", char[:]),
    ("numPlanets", int32),
    ("dim", int32),
    ("G", float64),
    ("pos", float64[:, :, :]),
    ("vel", float64[:, :, :])
]


#@jitclass(spec)
class solarsystem:
    def __init__(self, initPos, initVel, mass, N, tn, names):
        # Adding input
        self.mass = mass        # Array of mass, sequential
        self.N = int(N)         # Number of integration points
        self.tn = tn            # Simulation time [Yr]
        self.names = names      # Names of planets, for plotting
        # Some useful constants based on the input
        self.numPlanets = len(self.mass)
        self.h = self.tn / self.N
        self.dim = len(initPos[0][:])  # Number of dimensions, should be 2 or 3
        self.G = 4 * np.pi ** 2  # Gravitational Constant [AU^3 yr^-2 M_sun^-1]

        # Creating arrays
        self.pos = np.zeros([self.numPlanets, self.N, self.dim], dtype=np.float64)
        self.vel = np.zeros([self.numPlanets, self.N, self.dim], dtype=np.float64)
        # Setting initial conditions for all planets
        self.pos[:, 0] = initPos  # [AU]
        self.vel[:, 0] = initVel  # [AU/Yr]

        return

    def gravity(self, planetIndex, timeIndex):
        accel = np.zeros(self.dim)
        if planetIndex == 0:
            return accel
        else:
            pass

        for j in xrange(self.numPlanets):
            if j != planetIndex:  # Not interested in self-attraction
                relPos = self.pos[j, timeIndex] - self.pos[planetIndex, timeIndex]
                accel += relPos * self.G * self.mass[j] / np.linalg.norm(relPos) ** 3

        return accel

    def eulercromer(self, diffeq):
        """ 
        Solves N coupled differential equations using Euler-Cromer
        for a given differential equation
        """
        for i in xrange(self.N - 1):
            for n in xrange(self.numPlanets):
                acc = diffeq(n, i)
                self.vel[n][i + 1] = self.vel[n][i] + self.h * acc
                self.pos[n][i + 1] = self.pos[n][i] + self.h * self.vel[n][i + 1]

        return

    def velocityverlet(self, diffeq):
        for i in xrange(self.N - 1):
            for n in xrange(self.numPlanets):
                a1 = diffeq(n, i)       # Acceleration at t=t
                self.pos[n][i + 1] = self.pos[n][i] + self.vel[n][i] * self.h + 0.5 * a1 * self.h**2
                a2 = diffeq(n, i + 1)   # Acceleration at t=t+h
                self.vel[n][i + 1] = self.vel[n][i] + 0.5 * self.h * (a1 + a2)
        return

    def changeref(self, referencePlanet):
        """"
        Changes the frame of reference such that the planet with index=referencePlanet
        has pos=(0, 0, 0) and vel=(0, 0, 0) for all t.
        """
        for i in xrange(self.numPlanets):
            self.pos[i] = self.pos[i] - self.pos[referencePlanet]
            self.vel[i] = self.pos[i] - self.vel[referencePlanet]

        return

    def plot(self, fname):
        if self.dim == 3:
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            for i in range(self.numPlanets):
                # Data for a three-dimensional line
                xline = self.pos[i, :, 0]
                yline = self.pos[i, :, 1]
                zline = self.pos[i, :, 2]
                ax.plot3D(xline, yline, zline, label=self.names[i])
            ax.set_xlabel('x [AU]')
            ax.set_ylabel('y [AU]')
            ax.set_zlabel('z [AU]')
            plt.legend()
            plt.savefig("../figs/" + fname + ".pdf")
            plt.show()

        if self.dim == 2:
            for i in range(self.numPlanets):
                plt.plot(self.pos[i, :, 0], self.pos[i, :, 1], label=self.names[i])
            plt.xlabel("x [AU]")
            plt.ylabel("y [AU]")
            plt.legend()
            plt.savefig("../figs/" + fname + ".pdf")
            plt.show()

        return

    def get(self):
        return self.pos, self.vel


if __name__ == '__main__':
    m = np.array([3.00348959632E-6, 1])
    x0 = np.array([[1, 0], [0, 0]])
    v0 = np.array([[0, 2 * np.pi], [0, 0]])
    #planets = {"Earth": 399, "Sun": 10}
    #x01, v01, m1 = hori.fetch_data(jpl_id=planets)
    N = 1e5
    tn = 1
    names = ["Earth", "Sun"]

    esys = solarsystem(initPos=x0, initVel=v0, mass=m, N=N, tn=tn, names=names)
    esys.velocityverlet(esys.gravity)
    esys.plot("earthsun")
