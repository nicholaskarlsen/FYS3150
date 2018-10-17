# Uses astroquery to fetch JPL Horizons data to use as initial conditions & testing
import numpy as np
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


def fetch_data(jpl_id, dimensions, referenceFrame="500@0"):
    """ This function fetches position and velocity data from the Horizons system
    hosted by the Solar dynamics group at JPL, NASA using astroquery.
    Since i only need initial conditions, i fetch data using the default epoch
    kwarg in the Horizons() class, which yields the position and velocity vectors 
    at time = runtime """
    NumPlanets = len(dict)
    initPos = np.zeros([NumPlanets, dimensions])
    initVel = np.zeros([NumPlanets, dimensions])
    planetMass = np.zeros(NumPlanets)  # These have to be hardcoded

    for i, pname in enumerate(jpl_id):
        # Status update on a single, updating line
        print "\rFetching data from: https://ssd.jpl.nasa.gov/horizons_batch.cgi [%i/%i]" % (i, self.N_bodies),
        sys.stdout.flush()
        temp_obj = Horizons(id=jpl_id[pname], id_type='id',
                            location=referenceFrame)
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


if __name__ == '__main__':
    # Example showing how to use this script
    # List of planet names and their associated JPL Horizons ID
    planets = {"Sun": 10, "Mercury": 199, "Venus": 299, "Earth": 399, "Mars": 499,
               "Jupiter": 599, "Saturn": 699, "Uranus": 799, "Neptune": 899,
               "pluto": 999}
