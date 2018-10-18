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
    def __init__(self, initPos initVel, mass):
        self.numPlanets = len(jpl_planets)
        self.initPos = initPos
        self.initVel = initVel
        self.mass = mass

    def n_body_gravity(self, planet_index, time_index):
        """
        planet_index = Index of planet which is being solved
        time_index = Index of the time it is being solved at
        """
        accel = np.zeros(self.numPlanets - 1)  # N-1 planets acting on a particular one

        for i in xrange(self.numPlanets):
            if i != planet_index:  # Not interested in self-attraction
                rel_pos = (pos[planet_index] - pos)  # relative position vector
                a[i] = - self.G * mass[i] * rel_pos / np.linalg.norm(rel_pos)**3
        return sum(accel)

    @staticmethod
    @jit
    def Solver():
        return

    def test_function(self):
        test_planets = {"Sun": 10, "Earth": 399}

        planet = Horizons(id=planet_ID, id_type='id', location='500@0',
                          epochs={'start': '2018-10-11', 'stop': '2018-12-11',
                                  'step': '1y'})


if __name__ == '__main__':
    x0 = np.array([1, 0])  # (x, y) [AU]
    v0 = np.array([0, 4*np.pi])
