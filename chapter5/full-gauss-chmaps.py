from matplotlib import use, rc, cm, _cm
use('agg')
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/imaps/holly/home/ajl7/CoronaTemps/')
from temperature import TemperatureMap as tm
from sunpy.map import Map
from sunpy.time import parse_time as parse
from skimage import measure
from sunpy.net import hek
from os.path import join, expanduser

fs = 20
CThome = join(expanduser('~'), 'CoronaTemps')

#dates = ['2011-01-28', '2011-02-08', '2011-02-21'] # QS
#coords = [([200, 700], [-700, -200]), ([-400, 100], [-400, 100]), ([-500, 0], [-500, 0])] # QS
dates = ['2011-02-01', '2011-02-01', '2011-02-14'] # CH
coords = [([-500, 500], [0, 1000]), ([-400, 300], [-1150, -450]), ([-400, 400], [-1200, -400])] # CH
#dates = ['2011-01-22', '2011-01-20', '2011-02-01', '2011-02-19'] # AR
#coords = [([-200, 200], [300, 700]), ([-550, -250], [350, 650]),
#          ([-600, -200], [-450, -50]), ([-200, 200], [100, 500])] # AR

hfig, hax = plt.subplots(figsize=(16, 12))
for date, coord, thing, label, c in zip(dates, coords, ['a', 'b', ''], ['$CH_1$', '$CH_2$', '$CH_3$'], ['green', 'red', 'blue']):
    qsmap = tm(date, n_params=3,# verbose=True,
               data_dir=join(CThome, 'data'), maps_dir=CThome)
    emmap = Map(qsmap.emission_measure, qsmap.meta).submap(*coord)
    wmap = Map(qsmap.dem_width, qsmap.meta).submap(*coord)
    qsmap.save()
    qsmap = qsmap.submap(*coord)

    fig = plt.figure(figsize=(32, 24))
    qsmap.plot(vmin=5.8, vmax=6.5)
    #labels = ax1.get_xlabel(), ax1.get_ylabel()
    #ax1.set_xlabel(labels[0], fontsize=fs)
    #ax1.set_ylabel(labels[1], fontsize=fs)
    plt.title('Coronal hole {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        qsmap.date, np.nanmin(qsmap.data), np.nanmean(qsmap.data), np.nanmax(qsmap.data)))#,
    #    fontsize=24)
    plt.colorbar()
    
    plt.savefig('ch_{:%Y-%m-%dT%H%M}{}_full'.format(parse(qsmap.date), thing))
    plt.close()

    cmap = _cm.cubehelix(s=2.8, r=-0.7, h=1.4, gamma=1.0)
    cm.register_cmap(name='emhelix', data=cmap)
    emmap.cmap = cm.get_cmap('emhelix')
    fig = plt.figure(figsize=(32, 24))
    emmap.plot(vmin=np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)), 
               vmax=np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)))
    plt.title('Coronal hole $EM_{full}$')
    plt.colorbar()
    plt.savefig('chem_{:%Y-%m-%dT%H%M}{}_full'.format(parse(qsmap.date), thing))
    plt.close()
    print 'EM_full', label, np.nanmean(emmap.data)

    cmap = _cm.cubehelix(s=1.8, r=0.8, h=1.2, gamma=1.0)
    cm.register_cmap(name='whelix', data=cmap)
    wmap.cmap = cm.get_cmap('whelix')
    fig = plt.figure(figsize=(32, 24))
    wmap.plot(vmin=0.1, vmax=0.7)
    plt.title('Coronal hole DEM width')
    plt.colorbar()
    plt.savefig('chw_{:%Y-%m-%dT%H%M}{}_full'.format(parse(qsmap.date), thing))
    plt.close()
    print 'w', label, np.nanmean(wmap.data)

    plt.sca(hax)
    hist, bins = np.histogram(qsmap.data.flatten(), range=(5.6, 7.0), bins=71)
    hist = (hist / float(qsmap.size)) * 100
    plt.step(bins[:-1], hist, label=label, color=c)
    plt.xlabel('log(T)')
    plt.ylabel('% of image')

plt.legend()
plt.savefig('ch-histograms-full')
plt.close()
