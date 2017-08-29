from matplotlib import use, rc, cm, _cm
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='25.0')
import sys
from os.path import expanduser, join
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as tm
import matplotlib.pyplot as plt
from sunpy.map import Map
from sunpy.time import parse_time as pt
import numpy as np
from astropy import units as u

home = expanduser('~')
CThome = join(home, 'CoronaTemps')
datahome = join('/fastdata', 'sm1ajl', 'thesis')

"""
Plot temperatureand EM values using the full-Gaussian method.
"""

date = '2011-02-15'

# -----------
# Temperature
# -----------
tmap = tm(date, n_params=3, verbose=True,
          data_dir=join(datahome, 'data/2011/02/15'),
          maps_dir=join(datahome, 'maps/2011/02/15'))
tmap.save()

print tmap.min(), tmap.mean(), tmap.max()

fig = plt.figure(figsize=(32, 24))
tmap.plot(vmin=5.9, vmax=6.6)
plt.title('Full-Gaussian log(T) solution')
plt.colorbar()

plt.savefig(join(datahome, 'plots', 'chapter2', 't_3par_fd_2011-02-15'))
plt.close()
