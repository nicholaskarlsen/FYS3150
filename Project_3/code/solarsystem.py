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

# Default server does not work on my version of astroquery, so set it manually
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


class solarsystem:
    def __init__(self, initPos, initVel, mass, dt, tn):
        self.numPlanets = len(jpl_planets)
        self.initPos = initPos
        self.initVel = initVel
        self.mass = mas
        self.N = tn * dt
        self.pos = np.zeros(N)
        self.vel = np.zeros(N)
        self._G = 4 * np.pi**2    # Gravitational Constant [AU^3 yr^-2 M_sun^-1]

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
        accel = np.zeros(self.numPlanets - 1)  # N-1 planets acting on a particular one

        for i in xrange(self.numPlanets):
            if i != planet_index:  # Not interested in self-attraction
                rel_pos = (pos[planet_index] - pos[i])  # relative position vector
                accel[i] = accel_g(rel_pos, self.mass[i])
        return sum(accel)

    @staticmethod  # Needs to be static for JIT
    def eulercromer():
        for i in xrange(N - 1):
            for j in xrange():
            vel[i + 1] = vel[i] + h * acc
            pos[i + 1] = pos[i] + h * vel[i + 1]
        return

    @staticmethod  # Needs to be static for JIT
    def velocityverlet():
        for i in xrange(N - 1):
            a1 = diffeq(pos[i], mass=1)   # mass of sun = 1
            pos[i + 1] = pos[i] + vel[i] * h + 0.5 * a1 * h**2
            a2 = diffeq(pos[i + 1], mass=1)
            vel[i + 1] = vel[i] + 0.5 * h * (a1 + a2)
        return

    def plot(self):
        """
        Plots the data generated with the orbit method
        """
        for i in range(system.number_of_planets):
            pos = self.orbit(planet_index=i, time_stop=100, num_steps=1e6)[1]
            plt.plot(pos[:, 0], pos[:, 1], label="%d" % i)
        plt.title("Orbit of planets")
        plt.xlabel("$x$ [AU]")
        plt.ylabel("$y$ [AU]")
        plt.plot(0, 0, "r.")
        plt.legend()
        plt.annotate('Star', xy=(0, 0), xytext=(0, 0))
        # plt.savefig("orbitfig.png")
        plt.show()

        return

    def test_function(self):
        test_planets = {"Sun": 10, "Earth": 399}

        planet = Horizons(id=planet_ID, id_type='id', location='500@0',
                          epochs={'start': '2018-10-11', 'stop': '2018-12-11',
                                  'step': '1y'})


if __name__ == '__main__':
    from horizons import fetch_data

    sun_body_center = "500@10"  # JPL code for 'Sun (body center)' frame of reference

    planets = {"Earth": 399}
    x0, v0, m = fetch_data(jpl_id=planets, referenceFrame=sun_body_center)

    print x0, v0
