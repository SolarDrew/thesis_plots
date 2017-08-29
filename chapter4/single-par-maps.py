from matplotlib import use, rc
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='25.0')
import sys
import os
os.system('pwd')
from os.path import expanduser, join
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as tm
import matplotlib.pyplot as plt
from sunpy.map import Map
from sunpy.time import parse_time as pt
import numpy as np

home = expanduser('~')
CThome = join(home, 'CoronaTemps')
datahome = join('/fastdata', 'sm1ajl', 'thesis')
plotdir = join(datahome, 'plots', 'chapter4')

"""
Plot temperature and EM values using single-parameter methods.
"""

date = '2011-02-15'

# -----------
# Temperature
# -----------
tmap = tm(date, n_params=1, verbose=True,
          data_dir=join(datahome, 'data/{}/'.format(date).replace('-', '/')),
          maps_dir=join(datahome, 'maps/{}/'.format(date).replace('-', '/')))
tmap.save()

print tmap.min(), tmap.mean(), tmap.max()

fig = plt.figure(figsize=(32, 24))
tmap.plot(vmin=5.9, vmax=6.6)
plt.title('Single-parameter log(T) solution')
plt.colorbar()

plt.savefig(join(plotdir, 't_1par_fd_2011-02-15'))
plt.close()

"""fig = plt.figure(figsize=(32, 24))
tmap.plot()
plt.title('Single-parameter log(T) solution')
plt.colorbar()

plt.savefig('t_1par_fd_2011-02-15-struct')
plt.close()"""

# ----------------
# Emission Measure
# ----------------

# Individual wavelengths
alldata = np.zeros(shape=(4096, 4096))
threedata = np.zeros(shape=(4096, 4096))
for wl in ['94', '131', '171', '193', '211', '335']:
    emmap = tmap.calculate_em(wlen=wl)
    alldata += emmap.data
    if wl in ['171', '193', '211']:
        threedata += emmap.data

    print wl, np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)

    fig = plt.figure(figsize=(32, 24))
    """meanem = np.nanmean(emmap.data, dtype='float64')
    stdem = np.nanstd(emmap.data, dtype='float64')
    emmap.plot(vmin=meanem-(2.0*stdem),
               vmax=meanem+(2.0*stdem))"""
    emmap.plot(vmin=20.0, vmax=32.0)
    plt.title('Single-parameter log(EM) solution')
    plt.colorbar()

    plt.savefig(join(plotdir, 'em_1par-{}_fd_2011-02-15'.format(wl)))
    plt.close()

    emmap = emmap.submap([2000, 3024], [1200, 2224], units='pixels')
    print np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
    fig = plt.figure(figsize=(32, 24))
    emmap.plot(vmin=26.5, vmax=29.0, cmap='afmhot')
    plt.title('Single-parameter log(EM) solution')
    plt.colorbar()

    plt.savefig(join(plotdir, 'em_1par-{}_fd_2011-02-15_ar'.format(wl)))
    plt.close()


# Average over 171, 193 and 211
emmap.data = threedata / 3.0
print 'Three', np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
"""meanem = np.nanmean(emmap.data, dtype='float64')
stdem = np.nanstd(emmap.data, dtype='float64')
emmap.plot(vmin=meanem-(2.0*stdem),
           vmax=meanem+(2.0*stdem))"""
emmap.plot(vmin=20.0, vmax=32.0)
plt.title('Single-parameter log(EM) solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-three_fd_2011-02-15'.format(wl)))
plt.close()

emmap = emmap.submap([2000, 3024], [1200, 2224], units='pixels')
print np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
emmap.plot(vmin=26.5, vmax=29.0, cmap='afmhot')
plt.title('Single-parameter log(EM) solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-three_fd_2011-02-15_ar'.format(wl)))
plt.close()

# Average over all wavelengths
emmap.data = alldata / 6.0
print 'All', np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
"""meanem = np.nanmean(emmap.data, dtype='float64')
stdem = np.nanstd(emmap.data, dtype='float64')
emmap.plot(vmin=meanem-(2.0*stdem),
           vmax=meanem+(2.0*stdem))"""
emmap.plot(vmin=20.0, vmax=34.0)
plt.title('Single-parameter log(EM) solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-all_fd_2011-02-15'.format(wl)))
plt.close()

emmap = emmap.submap([2000, 3024], [1200, 2224], units='pixels')
print np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
emmap.plot(vmin=26.5, vmax=29.0, cmap='afmhot')
plt.title('Single-parameter log(EM) solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-all_fd_2011-02-15_ar'.format(wl)))
plt.close()

# ---------------
# Goodness-of-fit
# ---------------
gmap = Map(np.log10(tmap.goodness_of_fit), tmap.meta)
lowx, highx = (gmap.xrange[0] / gmap.scale['x'],
               gmap.xrange[1] / gmap.scale['x'])
lowy, highy = (gmap.yrange[0] / gmap.scale['y'],
               gmap.yrange[1] / gmap.scale['y'])
#x_grid, y_grid = np.mgrid[lowx:highx-1, lowy:highy-1]
x_grid, y_grid = np.mgrid[lowx:highx, lowy:highy]
r_grid = np.sqrt((x_grid ** 2.0) + (y_grid ** 2.0))
outer_rad = (gmap.rsun_arcseconds * 1.5) / gmap.scale['x']
print gmap.shape, x_grid.shape, y_grid.shape, r_grid.shape, outer_rad
gmap.data[r_grid > outer_rad] = None
print np.nanmin(gmap.data), np.nanmean(gmap.data, dtype='float64'), np.nanmax(gmap.data), np.nanstd(gmap.data, dtype='float64')

fig = plt.figure(figsize=(32, 24))
gmap.plot(cmap='cubehelix', vmin=-2.4, vmax=-0.6)
plt.title('Single-parameter solution goodness-of-fit')
plt.colorbar()

plt.savefig(join(plotdir, 'g_1par_fd_2011-02-15'))
plt.close()
