from matplotlib import use, rc, cm, _cm
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='25.0')
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import sys
sys.path.append('/home/sm1ajl/CoronaTemps/')
from temperature import TemperatureMap as tm
from sunpy.map import Map
from sunpy.time import parse_time as parse
from skimage import measure
from sunpy.net import hek
import os
from os.path import join, expanduser, exists

fs = 20
CThome = join(expanduser('~'), 'CoronaTemps')
datahome = join('/fastdata', 'sm1ajl', 'thesis', 'data')
mapshome = datahome.replace('/data', '/maps')
plotshome = join(datahome.replace('/data', '/plots'), 'chapter5', 'qs')
if not exists(plotshome):
    os.makedirs(plotshome)

qsmaps = []
dates = ['2011-01-28', '2011-02-08', '2011-02-21']
coords = [([200, 700], [-700, -200]), ([-400, 100], [-400, 100]), ([-500, 0], [-500, 0])]

for date, coord in zip(dates, coords):
    d = parse(date)
    fdmap = Map(join(datahome, '{0:%Y/%m/%d}/171/*_{0:%Y?%m?%d}?{0:%H?%M?%S}*fits').format(d))
    if isinstance(fdmap, list): fdmap = fdmap[0]
    fig, ax = plt.subplots(figsize=(16, 12))
    fdmap.plot()
    corner = coord[0][0], coord[1][0]
    x, y = coord[0][1] - coord[0][0], coord[1][1] - coord[1][0]
    rect = patches.Rectangle(corner, x, y, color='white', fill=False)
    ax.add_artist(rect)
    plt.savefig(join(plotshome, 'fulldisk_qs_{:%Y-%m-%dT%H%M}').format(d))
    plt.close()
    thismap = tm(date, verbose=True,
                 data_dir=join(datahome, '{:%Y/%m/%d}'.format(d)),
                 maps_dir=join(mapshome, '{:%Y/%m/%d}'.format(d)))
    thismap.save()
    qsmaps.append(thismap.submap(*coord))

hfig, hax = plt.subplots(figsize=(16, 12))
for qsmap, label, c in zip(qsmaps, ['$QS_1$', '$QS_2$', '$QS_3$'], ['green', 'red', 'blue']):
    fig = plt.figure(figsize=(32, 24))
    ax1 = fig.add_subplot(1, 1, 1)
    qsmap.plot(vmin=5.8, vmax=6.2)
    #labels = ax1.get_xlabel(), ax1.get_ylabel()
    #ax1.set_xlabel(labels[0], fontsize=fs)
    #ax1.set_ylabel(labels[1], fontsize=fs)
    plt.title('Quiet sun region {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        qsmap.date, np.nanmin(qsmap.data), np.nanmean(qsmap.data), np.nanmax(qsmap.data)))#,
    #    fontsize=24)
    plt.colorbar()
    
    plt.savefig(join(plotshome, 'qs_{:%Y-%m-%dT%H%M}').format(parse(qsmap.date)))
    plt.close()

    threedata = np.zeros(shape=qsmap.shape)
    alldata = np.zeros(shape=qsmap.shape)
    for wlen in ['94', '131', '171', '193', '211', '335']:
        fig = plt.figure(figsize=(32, 24))
        emmap = qsmap.calculate_em(wlen)
        alldata += emmap.data
        if wlen in ['171', '193', '211']: threedata += emmap.data
        #emmap.plot(vmin=np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
        #           vmax=np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)))
        emmap.plot(vmin=24.0, vmax=29.5)
        plt.title('Quiet sun $EM_{'+wlen+'}$')
        plt.colorbar()
        plt.savefig(join(plotshome, 'qs_{:%Y-%m-%dT%H%M}_em{}').format(parse(qsmap.date), wlen))
        plt.close()
    cmap = emmap.cmap

    fig = plt.figure(figsize=(32, 24))
    emmap = Map(threedata/3, qsmap.meta)
    #emmap.plot(vmin=np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
    #           vmax=np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)),
    #           cmap=cmap)
    emmap.plot(vmin=24.0, vmax=29.5, cmap=cmap)
    plt.title('Quiet sun $EM_{three}$')
    plt.colorbar()
    plt.savefig(join(plotshome, 'qs_{:%Y-%m-%dT%H%M}_emthree').format(parse(qsmap.date), wlen))
    plt.close()
    print 'EM_three', label, np.nanmean(emmap.data)

    fig = plt.figure(figsize=(32, 24))
    emmap = Map(alldata/6, qsmap.meta)
    #emmap.plot(vmin=np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
    #           vmax=np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)),
    #           cmap=cmap)
    emmap.plot(vmin=24.0, vmax=29.5, cmap=cmap)
    plt.title('Quiet sun $EM_{all}$')
    plt.colorbar()
    plt.savefig(join(plotshome, 'qs_{:%Y-%m-%dT%H%M}_emall').format(parse(qsmap.date), wlen))
    plt.close()

    plt.sca(hax)
    hist, bins = np.histogram(qsmap.data.flatten(), range=(5.6, 7.0), bins=71)
    hist = (hist / float(qsmap.size)) * 100
    plt.step(bins[:-1], hist, label=label, color=c)
    plt.xlabel('log(T)')
    plt.ylabel('% of image')

plt.legend()
plt.savefig(join(plotshome, 'qs-histograms'))
plt.close()
