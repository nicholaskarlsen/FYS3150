# Uses astroquery to fetch JPL Horizons data to use as initial conditions & testing
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

# List of planet names and their associated JPL Horizons ID
planets = {"Sun": 10, "Mercury": 199, "Venus": 299, "Earth": 399, "Mars": 499,
           "Jupiter": 599, "Saturn": 699, "Uranus": 799, "Neptune": 899,
           "pluto": 999}


ID = "Earth"
planet = Horizons(id=planets[ID], id_type='id', location='500@0',
                  epochs={'start': '2018-10-11', 'stop': '2018-12-11',
                  'step': '1y'})

print planet.vectors()[0]
