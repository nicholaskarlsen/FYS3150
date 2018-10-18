# Uses astroquery to fetch JPL Horizons data to use as initial conditions & testing
import numpy as np
import sys
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


def fetch_data(jpl_id, referenceFrame="500@0"):
    """ This function fetches position and velocity data from the Horizons system
    hosted by the Solar dynamics group at JPL, NASA using astroquery.
    Since i only need initial conditions, i fetch data using the default epoch
    kwarg in the Horizons() class, which yields the position and velocity vectors 
    at time = runtime """
    NumPlanets = len(jpl_id)
    initPos = np.zeros([NumPlanets, 3])
    initVel = np.zeros([NumPlanets, 3])
    planetMass = np.zeros(NumPlanets)

    # Astroquery doesn't seemt to have masses, so i hardcode these
    # With values from the problem text [Kg]
    MASS_TABLE = {"Sun": 2e30, "Mercury": 3.3e23, "Venus": 4.9e24, "Earth": 6e24, "Mars": 6.6e23,
                  "Jupiter": 1.9e27, "Saturn": 5.5e26, "Uranus": 8.8e25, "Neptune": 1.03e26,
                  "Pluto": 1.31e22}

    for i, pname in enumerate(jpl_id):
        # Status update on a single, updating line
        print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [%i/%i]" % (i, NumPlanets),
        sys.stdout.flush()
        temp_obj = Horizons(id=jpl_id[pname], id_type='id',
                            location=referenceFrame)
        """Fetches position and velocity data in order (x, y, z)
        Note: Print method in vectors() doesnt seem to play nice with list
        comprehensions Hence the ugly (but stable and functioning)
        implemetation here."""
        initPos[i] = (temp_obj.vectors()["x"], temp_obj.vectors()["y"],
                      temp_obj.vectors()["z"])  # [AU]
        initVel[i] = (temp_obj.vectors()["vx"], temp_obj.vectors()["vy"],
                      temp_obj.vectors()["vz"])  # [AU/day]
        planetMass[i] = MASS_TABLE[pname]   # Fetches the mass from the hardcoded table
    print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [COMPLETE]"

    return initPos, initVel, planetMass


if __name__ == '__main__':
    # Example showing how to use this script
    # List of planet names and their associated JPL Horizons ID
    planets = {"Sun": 10, "Mercury": 199, "Venus": 299, "Earth": 399, "Mars": 499,
               "Jupiter": 599, "Saturn": 699, "Uranus": 799, "Neptune": 899,
               "Pluto": 999}

    x0, v0, m = fetch_data(planets, 3)
