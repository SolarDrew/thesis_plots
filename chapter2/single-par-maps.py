from matplotlib import use, rc
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
import sys
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
plotdir = join(datahome, 'plots', 'chapter2')

"""
Plot temperature and EM values using single-parameter methods.
"""

date = '2011-02-15'

tmap = tm(date, n_params=1, verbose=True,
          data_dir=join(datahome, 'data'),
          maps_dir=join(datahome, 'maps'))
tmap.save()

# Individual wavelengths
alldata = np.zeros(shape=(6, 4096, 4096))
threedata = np.zeros(shape=(3, 4096, 4096))
channels = ['94', '131', '171', '193', '211', '335']
w = 0
for wl in channels:
    emmap = tmap.calculate_em(wlen=wl)
    alldata[channels.index(wl)] = emmap.data
    if wl in ['171', '193', '211']:
        threedata[w] = emmap.data
        w += 1

    print wl, np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)

    fig = plt.figure(figsize=(32, 24))
    emmap.plot(vmin=20.0, vmax=34.0)
    plt.title('Single-parameter log(EM) solution')
    plt.colorbar()

    plt.savefig(join(plotdir, 'em_1par-{}_fd_2011-02-15'.format(wl)))
    plt.close()

# Average over 171, 193 and 211
emmap.data = np.nanmean(threedata, axis=0)# / 3.0
print 'Three', np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
emmap.plot(vmin=20.0, vmax=34.0)
plt.title('Single-parameter log(EM) solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-three_fd_2011-02-15'.format(wl)))
plt.close()

# Standard deviation of 171, 193 and 211
emmap.data = np.nanstd(threedata, axis=0)# / 3.0
print 'Three std. dev.', np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
d = emmap.data
emmap.plot(vmin=np.nanmean(d, dtype='float64')-(2*np.nanstd(d, dtype='float64')),
           vmax=np.nanmean(d, dtype='float64')+(2*np.nanstd(d, dtype='float64')))
plt.title('Single-parameter log(EM) solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-three-std_fd_2011-02-15'.format(wl)))
plt.close()

# Average over all wavelengths
emmap.data = np.nanmean(alldata, axis=0)# / 6.0
print 'All', np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
emmap.plot(vmin=20.0, vmax=34.0)
plt.title('Standard deviation of EM solutions')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-all_fd_2011-02-15'.format(wl)))
plt.close()

# Average over all wavelengths
emmap.data = np.nanstd(alldata, axis=0)# / 6.0
print 'All std. dev.', np.nanmin(emmap.data), np.nanmean(emmap.data, dtype='float64'), np.nanmax(emmap.data)
fig = plt.figure(figsize=(32, 24))
d = emmap.data
emmap.plot(vmin=np.nanmean(d, dtype='float64')-(2*np.nanstd(d, dtype='float64')),
           vmax=np.nanmean(d, dtype='float64')+(2*np.nanstd(d, dtype='float64')))
plt.title('Standard deviation of EM solution')
plt.colorbar()

plt.savefig(join(plotdir, 'em_1par-all-std_fd_2011-02-15'.format(wl)))
plt.close()
