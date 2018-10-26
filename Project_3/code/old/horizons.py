# Uses astroquery to fetch JPL Horizons data to use as initial conditions & testing
import numpy as np
import sys
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


def fetch_data(jpl_id, referenceFrame="500@0"):
    """
    This function fetches position and velocity data from the Horizons system
    hosted by the Solar dynamics group at JPL, NASA using astroquery.
    Since i only need initial conditions, i fetch data using the default epoch
    kwarg in the Horizons() class, which yields the position and velocity vectors
    at time = runtime
    """
    NumPlanets = len(jpl_id)
    initPos = np.zeros([NumPlanets, 3])
    initVel = np.zeros([NumPlanets, 3])
    planetMass = np.zeros(NumPlanets)

    # Astroquery doesn't seemt to have masses, so i hardcode these
    # Values collected from https://en.wikipedia.org/wiki/Planetary_mass
    # in units [Solar mass]
    MASS_TABLE = {"Sun": 1, "Mercury": 0.16601E-6, "Venus": 2.4478383E-6,
                  "Earth": 3.00348959632E-6, "Mars": 0.3227151E-6, "Jupiter": 954.79194E-6,
                  "Saturn": 285.8860E-6, "Uranus": 43.66244E-6, "Neptune": 51.51389E-6,
                  "Pluto": 0.007396E-6}


    for i, pname in enumerate(jpl_id):
        print "Here", i, pname
        # Status update on a single, updating line
        print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [%i/%i]" % (i, NumPlanets),
        sys.stdout.flush()
        temp_obj = Horizons(id=jpl_id[pname], id_type='id',
                            location=referenceFrame)
        """Fetches position and velocity data in order (x, y, z)
        Note: Print method in vectors() doesnt seem to play nice with list
        comprehensions Hence the ugly (but stable and functioning)
        implemetation here."""
        initPos[i] = [temp_obj.vectors()["x"], temp_obj.vectors()["y"],
                      temp_obj.vectors()["z"]]  # [AU]
        initVel[i] = [temp_obj.vectors()["vx"], temp_obj.vectors()["vy"],
                      temp_obj.vectors()["vz"]]  # [AU/day]
        initVel = initVel * (365.25)  # Convert to units [AU/yr]
        planetMass[i] = MASS_TABLE[pname]   # Fetches the mass from the hardcoded table
    print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [COMPLETE]"

    return initPos, initVel, planetMass


if __name__ == '__main__':
    # Example showing how to use this script
    # List of planet names and their associated JPL Horizons ID
    planets = {"Sun": 10, "Mercury": 199, "Venus": 299, "Earth": 399, "Mars": 499,
               "Jupiter": 599, "Saturn": 699, "Uranus": 799, "Neptune": 899,
               "Pluto": 999}

    x0, v0, m = fetch_data(planets)
