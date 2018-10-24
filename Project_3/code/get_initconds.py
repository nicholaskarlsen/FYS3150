# Uses astroquery to fetch JPL Horizons data to use as initial conditions & testing
import numpy as np
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'


def get_data(id_list, referenceFrame="500@10"):   # defaut @ sun center of mass

    initpos = np.zeros([len(id_list), 3], dtype=np.float64)
    initvel = np.zeros([len(id_list), 3], dtype=np.float64)

    for i in xrange(len(id_list)):
        obj = Horizons(id=id_list[i], id_type='id',
                       epochs={'start': '2018-10-16', 'stop': '2018-11-15', 'step': '1d'},
                       location=referenceFrame)
        vec = obj.vectors()
        x = vec['x'][0]
        y = vec['y'][0]
        z = vec['z'][0]
        vx = vec['vx'][0]
        vy = vec['vy'][0]
        vz = vec['vz'][0]

        initpos[i] = [x, y, z]      # [AU]
        initvel[i] = [vx, vy, vz]   # [AU/Day]

    initvel *= 365.25   # [AU/Day] -> [AU/Year]

    return initpos, initvel


if __name__ == '__main__':
    '*** How to use ****'
    # Make list of planet ids (see JPL Horizons website)
    lst = [10, 399, 599]
    # Call get_data with lst parameter and get initial conditions
    r0, v0 = get_data(lst)
