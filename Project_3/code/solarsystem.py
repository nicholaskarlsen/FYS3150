#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Nicholas Karlsen
# Python version 2.7.15
"""
The class fetches its data from the JPL Horizons database during runtime using astroquery,
therefore an active internet connection is required, however, the script can easily
be modified to fetch and save and offline copy if so desired.
"""
from __future__ import division  # Nobody expects the integer division
import numpy as np
from astroquery.jplhorizons import Horizons
import solver
import time
import sys
from astroquery.jplhorizons import conf

conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'  # Default does not work on my version of astroquery


class solarsystem:
    def __init__(self, jpl_planets, Dim=3):
        self.jpl_planets = jpl_planets
        self.N_bodies = len(jpl_planets)
        self.Dim = Dim

    def fetch_data(self):
        """Fetches data from JPL Horizons"""
        initPos = np.zeros([self.N_bodies, self.Dim])
        initVel = np.zeros([self.N_bodies, self.Dim])
        planetMass = np.zeros(self.N_bodies)

        print "Fetching data"
        for i, pname in enumerate(self.jpl_planets):
            temp_obj = Horizons(id=self.jpl_planets[pname], id_type='id', location='500@0',
                                epochs={'start': '2014-10-11', 'stop': '2016-12-11', 'step': '1y'})
            # Fetches position and velocity data in order (x, y, z)
            initPos[i] = (temp_obj.vectors()[0][3], temp_obj.vectors()[0][4], temp_obj.vectors()[0][5])  # [AU]
            initVel[i] = (temp_obj.vectors()[0][6], temp_obj.vectors()[0][7], temp_obj.vectors()[0][8])  # [AU/day]
        print "Done"
        return

    def test_function(self):
        test_planets = {"Sun": 10, "Earth": 399}

        planet = Horizons(id=planet_ID, id_type='id', location='500@0',
                          epochs={'start': '2018-10-11', 'stop': '2018-12-11',
                                  'step': '1y'})


if __name__ == '__main__':
    planets = {"Sun": 10, "Earth": 399}
    test = solarsystem(planets)
    test.fetch_data()
