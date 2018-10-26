#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15
"""
The class fetches its data from the JPL Horizons database during runtime using 
astroquery, therefore an active internet connection is required, however, the script 
can easily be modified to fetch and save and offline copy if so desired.

The class takes an arguement, that is a dictionary formatted such that each entry is
{'name string':id int}
Name string isn't strict, and can be of your choosing, but the planet name makes sense.
id string needs to correspond EXACTLY to the id number which is used in the jpl horizons system
"""
from __future__ import division  # Nobody expects the integer division
import numpy as np
from numba import jit
from astroquery.jplhorizons import Horizons, conf
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from numba import jitclass          # import the decorator
from numba import int32, float32    # import the types

# Default server does not work on my version of astroquery, so set it manually
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


class solarsystem:
    def __init__(self, initPos, initVel, mass, N, tn):
        """
        The class is initialized by feeding it arrays containing relevant info
        pertaining to the system at hand.
        initPos : [planet, time, dimension]
        initVel : [planet, time, dimension]
        mass : [planet]
        N : Number of integration points
        """

        # First check that the shape of input arrays make sense
        if np.size(initPos) != np.size(initVel):
            raise ValueError("shape of initPos != shape of initVel")

        if len(initPos) != len(mass):
            raise ValueError("Length of initPos/initVel != length of mass array.")

        # Define some useful numbers
        self._G = 4 * np.pi**2    # Gravitational Constant [AU^3 yr^-2 SolarMass^-1]
        self.N = int(N)          # Number of integration points
        self.tn = tn        # Simulation time / time at N. [Yr]
        self.dt = tn / N    # Timestep
        self.numPlanets = len(initPos)

        self.mass = mass    # Mass of planets, array. [SolarMass]
        # Store position & velocity data for ALL planets in 2 arrays
        self.pos = np.zeros([self.numPlanets, self.N, 3])  # [AU]
        self.vel = np.zeros([self.numPlanets, self.N, 3])  # [AU/Yr]
        # Setting initial conditions for all planets
        self.pos[:, 0] = initPos
        self.vel[:, 0] = initVel
        self.dim = 3

    def accel_g(self, position_vector, mass):
        '''
        Returns the gravitational acc from the star for a given vector (i.e works for 1,
        2 and 3-Dim)
        '''
        return - pos * self.G * mass / np.linalg.norm(pos) ** 3

    def n_body_gravity(self, planet_index, time_index):
        """
        planet_index = Index of planet which is being solved
        time_index = Index of the time it is being solved at
        """
        accel = np.zeros(self.dim)

        for i in xrange(self.numPlanets):
            if i != planet_index:  # Not interested in self-attraction
                # compute relative position vector, then acceleration due to that body
                rel_pos = (self.pos[planet_index, time_index] - self.pos[i, time_index])
                # print rel_pos, i
                accel += rel_pos * self._G * self.mass[i] / np.linalg.norm(rel_pos) ** 3
        return accel  # Return sum of acceleration due to all bodies.

    def eulercromer(self, diffeq):
        """Static methods doesn't have access to self variables, so these need to be"""
        for i in xrange(self.N - 1):
            for j in xrange(self.numPlanets):
                acc = diffeq(j, i)
                self.vel[j, i + 1] = self.vel[j, i] + self.dt * acc
                self.pos[j, i + 1] = self.pos[j, i] + self.dt * self.vel[j, i + 1]
        return

    def velocityverlet(self, diffeq):
        for i in xrange(self.N - 1):
            for j in xrange(self.numPlanets):
                a1 = diffeq(j, i)   # mass of sun = 1
                pos[i + 1] = pos[i] + vel[i] * h + 0.5 * a1 * h**2
                a2 = diffeq(j, i + 1)
                vel[i + 1] = vel[i] + 0.5 * h * (a1 + a2)
        return

    def plot(self):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for i in range(self.numPlanets):
            # Data for a three-dimensional line
            xline = self.pos[i][:, 0]
            yline = self.pos[i][:, 1]
            zline = self.pos[i][:, 2]
            ax.plot3D(xline, yline, zline)
        ax.set_xlabel('x [AU]')
        ax.set_ylabel('y [AU]')
        ax.set_zlabel('z [AU]')
        plt.show()

        return

    def get(self):
        return self.pos, self.vel

    def test_function(self):
        test_planets = {"Sun": 10, "Earth": 399}

        planet = Horizons(id=planet_ID, id_type='id', location='500@0',
                          epochs={'start': '2018-10-11', 'stop': '2018-12-11',
                                  'step': '1y'})

if __name__ == '__main__':
    from horizons import fetch_data

    planets = {"Earth": 399, "Sun": 10, "Jupiter": 599}
    x0, v0, m = fetch_data(jpl_id=planets)

    mysys = solarsystem(x0, v0, m, 1e7, 1)

    mysys.eulercromer(mysys.n_body_gravity)
    mysys.plot()
