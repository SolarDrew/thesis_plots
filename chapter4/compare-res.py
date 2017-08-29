import numpy as np
from matplotlib import use, rc
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='25.0')
import matplotlib.pyplot as plt
from matplotlib import patches
import sunpy
from sunpy.map import Map
from os import path, makedirs
import sys
CThome = path.join(path.expanduser('~'), 'CoronaTemps')
sys.path.append(CThome)
from temperature import TemperatureMap as tm
from utils import gaussian, load_temp_responses

datahome = path.join('/fastdata', 'sm1ajl', 'thesis', 'data')
mapshome = datahome.replace('/data', '/maps')
plotshome = path.join(mapshome.replace('maps', 'plots'), 'chapter4')

date = sunpy.time.parse_time('2011-02-15')

single = tm(date,
            data_dir=path.join(datahome, '{:%Y/%m/%d}'.format(date)),
            maps_dir=path.join(mapshome, '{:%Y/%m/%d}'.format(date))).submap(
    [100, 300],
    [-300, -100])
sEM = single.calculate_em('171')

fullG = tm(date, n_params=3,
           data_dir=path.join(datahome, '{:%Y/%m/%d}'.format(date)),
           maps_dir=path.join(mapshome, '{:%Y/%m/%d}'.format(date))).submap(
    [100, 300],
    [-300, -100])

c = 175, -240

x, y = int(single.data_to_pixel(c[0], 'x')), int(single.data_to_pixel(c[1], 'y'))
print x, y

fig, ax = plt.subplots(1, 2, figsize=(28, 36))
plt.sca(ax[0])
single.plot()
ax[0].plot([c[0]], [c[1]], 'o', color='blue')
plt.colorbar(orientation='horizontal')
plt.sca(ax[1])
fullG.plot()
ax[1].plot([c[0]], [c[1]], 'o', color='blue')
plt.colorbar(orientation='horizontal')
plt.savefig(path.join(plotshome, 'bothmaps'))
plt.close()

logt = np.arange(5.6, 7.01, 0.01)
sgauss = gaussian(logt, single.data[x, y], 0.1, sEM.data[x, y])
fgauss = gaussian(logt, fullG.data[x, y], fullG.dem_width[x, y], fullG.emission_measure[x, y])
fig = plt.figure(figsize=(21, 18))
plt.plot(logt, sgauss)
plt.plot(logt, fgauss)
plt.savefig(path.join(plotshome, 'bothdems'))
plt.close()
