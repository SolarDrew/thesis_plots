from matplotlib import use, rc
use('agg')
rc('savefig', bbox='tight', pad_inches=0.5)
import sys
from os.path import expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as tm
from utils import load_temp_responses
import matplotlib.pyplot as plt
from sunpy.map import Map
from sunpy.time import parse_time as pt
import numpy as np
import glob
from calcemiss import reconstruct

""" 
Feed calculated temperatures back through response functions to create synthetic AIA images
"""

date = '2011-02-15'

# -----------------------------
# Single-parameter temperatures
# -----------------------------
tmap = tm(date, n_params=1,# verbose=True,
          data_dir=expanduser('~/CoronaTemps/data'),
          maps_dir=expanduser('~/CoronaTemps'))

logt = np.arange(0, 15.05, 0.05)
delta_t = logt[1] - logt[0]

resp = load_temp_responses()
wlens = ['94']#, '131', '171', '193', '211', '335']

for wl, thisresp in zip(wlens, resp):
    emmap = tmap.calculate_em(wl)
    if wl == '94': wl = '094'
    origfname = expanduser('~/CoronaTemps/data/{}/{}/*'.format(date.replace('-', '/'), wl))
    origfname = glob.glob(origfname)
    origmap = Map(origfname[-1])

    #thisresp /= thisresp.max()
    #emiss /= emiss.max()

    print thisresp.shape, thisresp.dtype
    emiss = reconstruct(tmap.data, thisresp, logt, tmap.shape[0], tmap.shape[1])
    thismap = Map(emiss, origmap.meta)

    fig = plt.figure(figsize=(64, 24))
    ax = fig.add_subplot(1, 2, 1)
    origmap.plot()
    plt.title('Original {} emission'.format(wl))
    plt.colorbar()
    ax = fig.add_subplot(1, 2, 2)
    thismap.plot()
    plt.title('{} emission inferred from single-parameter DEMs'.format(wl))
    plt.colorbar()

    plt.savefig('{}_1par_fd_2011-02-15'.format(wl))
    plt.close()
