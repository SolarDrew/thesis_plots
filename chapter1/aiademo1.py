from matplotlib import use, rc
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='25.0')
import matplotlib.pyplot as plt
import sunpy.map
from os import path

fig, axes = plt.subplots(3, 2, figsize=(32, 48))
axes = axes.flatten()

for i, wl in enumerate('94 131 171 193 211 335'.split()):
    if wl == '94':
        wl2 = '094'
    else:
        wl2 = wl
    fname = '/fastdata/sm1ajl/thesis/data/2011/02/15/{0}/*{1}*00_00*'.format(wl2, wl)
    print fname
    thismap = sunpy.map.Map(fname)
    if isinstance(thismap, list): thismap = thismap[0]
    plt.sca(axes[i])
    thismap.plot()

plt.savefig('/fastdata/sm1ajl/thesis/plots/chapter1/aiademo1')
plt.close()
