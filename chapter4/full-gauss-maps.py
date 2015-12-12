from matplotlib import use, rc, cm, _cm
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
import sys
from os.path import expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as tm
import matplotlib.pyplot as plt
from sunpy.map import Map
from sunpy.time import parse_time as pt
import numpy as np
from astropy import units as u

"""
Plot temperatureand EM values using the full-Gaussian method.
"""

date = '2011-02-15'

# -----------
# Temperature
# -----------
tmap = tm(date, n_params=3, verbose=True,
          data_dir=expanduser('~/CoronaTemps/data'),
          maps_dir=expanduser('~/CoronaTemps/maps'))
tmap.save()

print tmap.min(), tmap.mean(), tmap.max()

fig = plt.figure(figsize=(32, 24))
tmap.plot(vmin=5.9, vmax=6.6)
plt.title('Full-Gaussian log(T) solution')
plt.colorbar()

plt.savefig('t_3par_fd_2011-02-15')
plt.close()

fig = plt.figure(figsize=(32, 24))
tmap.plot()
plt.title('Full-Gaussian log(T) solution')
plt.colorbar()

plt.savefig('t_3par_fd_2011-02-15-struct')
plt.close()

# ----------------
# Emission Measure
# ----------------
emmap = Map(tmap.emission_measure, tmap.meta)

print np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)

cmap = _cm.cubehelix(s=2.8, r=-0.7, h=1.4, gamma=1.0)
cm.register_cmap(name='emhelix', data=cmap)
emmap.cmap = cm.get_cmap('emhelix')
fig = plt.figure(figsize=(32, 24))
#meanem = np.nanmean(emmap.data, dtype='float64')
#stdem = np.nanstd(emmap.data, dtype='float64')
#emmap.plot(vmin=meanem-(2.0*stdem),
#           vmax=meanem+(2.0*stdem))
emmap.plot(vmin=20.0, vmax=34.0)
plt.title('Full-Gaussian log(EM) solution')
plt.colorbar()

plt.savefig('em_3par_fd_2011-02-15')
plt.close()

emmap = emmap.submap([2000, 3024], [1200, 2224], units='pixels')
print np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
emmap.plot(vmin=26.5, vmax=29.0, cmap='afmhot')
plt.title('Full-Gaussian log(EM) solution')
plt.colorbar()

plt.savefig('em_3par_fd_2011-02-15_ar')
plt.close()

# -----
# Width
# -----
wmap = Map(tmap.dem_width, tmap.meta)
lowx, highx = (wmap.xrange[0] / wmap.scale['x'],
               wmap.xrange[1] / wmap.scale['x'])
lowy, highy = (wmap.yrange[0] / wmap.scale['y'],
               wmap.yrange[1] / wmap.scale['y'])
x_grid, y_grid = np.mgrid[lowx:highx-1, lowy:highy-1]
r_grid = np.sqrt((x_grid ** 2.0) + (y_grid ** 2.0))
outer_rad = (wmap.rsun_arcseconds * 1.5) / wmap.scale['x']
wmap.data[r_grid > outer_rad] = None
print np.nanmin(wmap.data), np.nanmean(wmap.data, dtype='float64'), np.nanmax(wmap.data), np.nanstd(wmap.data, dtype='float64')

cmap = _cm.cubehelix(s=1.8, r=0.8, h=1.2, gamma=1.0)
cm.register_cmap(name='whelix', data=cmap)
wmap.cmap = cm.get_cmap('whelix')
fig = plt.figure(figsize=(32, 24))
wmap.plot(vmin=0.0, vmax=0.7)
plt.title('Full-Gaussian solution goodness-of-fit')
plt.colorbar()

plt.savefig('w_3par_fd_2011-02-15')
plt.close()

# ---------------
# Goodness-of-fit
# ---------------
gmap = Map(np.log10(tmap.goodness_of_fit), tmap.meta)
lowx, highx = (gmap.xrange[0] / gmap.scale['x'],
               gmap.xrange[1] / gmap.scale['x'])
lowy, highy = (gmap.yrange[0] / gmap.scale['y'],
               gmap.yrange[1] / gmap.scale['y'])
x_grid, y_grid = np.mgrid[lowx:highx, lowy:highy]
r_grid = np.sqrt((x_grid ** 2.0) + (y_grid ** 2.0))
outer_rad = (gmap.rsun_arcseconds * 1.5) / gmap.scale['x']
gmap.data[r_grid > outer_rad] = None
print np.nanmin(gmap.data), np.nanmean(gmap.data, dtype='float64'), np.nanmax(gmap.data), np.nanstd(gmap.data, dtype='float64')

fig = plt.figure(figsize=(32, 24))
vmin = np.nanmean(gmap.data, dtype='float64') - (2.0 * np.nanstd(gmap.data, dtype='float64'))
vmax = np.nanmean(gmap.data, dtype='float64') + (2.0 * np.nanstd(gmap.data, dtype='float64'))
gmap.plot(cmap='cubehelix', vmin=vmin, vmax=vmax)
plt.title('Full-Gaussian solution goodness-of-fit')
plt.colorbar()

plt.savefig('g_3par_fd_2011-02-15')
plt.close()
