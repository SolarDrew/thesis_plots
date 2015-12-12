from matplotlib import use, rc
use('agg')
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import sunpy.map
from os import path

fig, axes = plt.subplots(2, 2, figsize=(32, 32))
axes = axes.flatten()

for i, wl in enumerate('304 1600 1700 4500'.split()):
#  try:
    fname = '/imaps/sspfs/archive/sdo/aia/fulldisk/data/2011/02/15/{0}/*{0}*00_00*'.format(wl)
    thismap = sunpy.map.Map(fname)
    if isinstance(thismap, list): thismap = thismap[0]
    print fname, thismap.exposure_time
    plt.sca(axes[i])
    thismap.plot()
#  except:
#    print 'Balls!', wl

plt.savefig('aiademo2')
plt.close()
