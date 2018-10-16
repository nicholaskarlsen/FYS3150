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
from astroquery.jplhorizons import Horizons, conf
import sys

# Default server does not work on my version of astroquery, so set it manually
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


class solarsystem:
    def __init__(self, jpl_planets, Dim=3):
        """
        Fetches data from JPL Horizons and places it in arrays let epoch
        arguement in Horizons() be default, i.e fetch data with current
        poisitions of planets during runtime -> init conditions will change
        every time the script is ran (unless you wait for the stars to align
        , again)
        """
        self.jpl_planets = jpl_planets
        self.N_bodies = len(jpl_planets)
        self.Dim = Dim
        initPos = np.zeros([self.N_bodies, self.Dim])
        initVel = np.zeros([self.N_bodies, self.Dim])
        planetMass = np.zeros(self.N_bodies)
        for i, pname in enumerate(self.jpl_planets):
            # Status update on a single, updating line
            print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [%i/%i]" % (i, self.N_bodies),
            sys.stdout.flush()
            temp_obj = Horizons(id=self.jpl_planets[pname], id_type='id',
                                location='500@0')
            """Fetches position and velocity data in order (x, y, z)
            Note: Print method in vectors() doesnt seem to play nice with list
            comprehensions Hence the ugly (but stable and functioning)
            implemetation here."""
            initPos[i] = (temp_obj.vectors()[0][3], temp_obj.vectors()[
                          0][4], temp_obj.vectors()[0][5])  # [AU]
            initVel[i] = (temp_obj.vectors()[0][6], temp_obj.vectors()[
                          0][7], temp_obj.vectors()[0][8])  # [AU/day]
        print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [COMPLETE]"
        return

    def gravityforce(self):
        """ Calculates the sum of all gravitational forces acting a particular
        body and returns the resultant Force vector """

        resultantForce =

        # pos_j - pos_i, i = 1,..., N,  j!=i
        for i, pos in enumerate(pos):
            if i != j:
                rel_pos = pos_j - pos[i]

        return

    def test_function(self):
        test_planets = {"Sun": 10, "Earth": 399}

        planet = Horizons(id=planet_ID, id_type='id', location='500@0',
                          epochs={'start': '2018-10-11', 'stop': '2018-12-11',
                                  'step': '1y'})


if __name__ == '__main__':
    planets = {"Sun": 10, "Earth": 399}
    test = solarsystem(planets)
