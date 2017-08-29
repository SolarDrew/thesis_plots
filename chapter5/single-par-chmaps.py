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
plotshome = join(datahome.replace('/data', '/plots'), 'chapter5', 'ch')
if not exists(plotshome):
    os.makedirs(plotshome)


dates = ['2011-02-01', '2011-02-01', '2011-02-14']
coords = [([-500, 500], [0, 1000]), ([-400, 300], [-1150, -450]), ([-400, 400], [-1200, -400])]
chmaps = []

for date, coord in zip(dates, coords):
    d = parse(date)
    if coord[0][1] != 300:
        fdmap = Map(join(datahome, '{0:%Y/%m/%d}/171/*_{0:%Y?%m?%d}?{0:%H?%M?%S}*fits').format(parse(date)))
        if isinstance(fdmap, list): fdmap = fdmap[0]
        fig, ax = plt.subplots(figsize=(16, 12))
        fdmap.plot()
    corner = coord[0][0], coord[1][0]
    x, y = coord[0][1] - coord[0][0], coord[1][1] - coord[1][0]
    rect = patches.Rectangle(corner, x, y, color='white', fill=False)
    ax.add_artist(rect)
    if coord[0][0] != -500:
        plt.savefig(join(plotshome, 'fulldisk_ch_{:%Y-%m-%dT%H%M}').format(parse(date)))
        plt.close()
    thismap = tm(date, verbose=True,
                 data_dir=join(datahome, '{:%Y/%m/%d}'.format(d)),
                 maps_dir=join(mapshome, '{:%Y/%m/%d}'.format(d)))
    thismap.save()
    chmaps.append(thismap.submap(*coord))

hfig, hax = plt.subplots(figsize=(16, 12))
for chmap, thing, label, c in zip(chmaps, ['a', 'b', ''], ['$CH_1$', '$CH_2$', '$CH_3$'], ['green', 'red', 'blue']):
    fig = plt.figure(figsize=(32, 24))
    ax1 = fig.add_subplot(1, 1, 1)
    chmap.plot()#vmin=5.8, vmax=6.5)
    #labels = ax1.get_xlabel(), ax1.get_ylabel()
    #ax1.set_xlabel(labels[0], fontsize=fs)
    #ax1.set_ylabel(labels[1], fontsize=fs)
    #chmap.draw_limb(color='black')
    #plt.axvline(x=0, ymin=0, ymax=0.375, color='blue')
    plt.title('Coronal hole {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        chmap.date, np.nanmin(chmap.data), np.nanmean(chmap.data), np.nanmax(chmap.data)))#,
    #    fontsize=24)
    plt.colorbar()
    
    plt.savefig(join(plotshome, 'ch_{:%Y-%m-%dT%H%M}{}').format(parse(chmap.date), thing))
    plt.close()

    vmin = 24.0
    vmax = 29.5
    threedata = np.zeros(shape=chmap.shape)
    alldata = np.zeros(shape=chmap.shape)
    for wlen in ['94', '131', '171', '193', '211', '335']:
        fig = plt.figure(figsize=(32, 24))
        emmap = chmap.calculate_em(wlen)
        alldata += emmap.data
        if wlen in ['171', '193', '211']: threedata += emmap.data
        #emmap.plot(vmin=np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
        #           vmax=np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)))
        emmap.plot(vmin=vmin, vmax=vmax)
        plt.title('Coronal hole $EM_{'+wlen+'}$')
        plt.colorbar()
        plt.savefig(join(plotshome, 'ch_{:%Y-%m-%dT%H%M}{}_em{}').format(parse(chmap.date), thing, wlen))
        plt.close()
    cmap = emmap.cmap

    fig = plt.figure(figsize=(32, 24))
    emmap = Map(threedata/3, chmap.meta)
    emmap.plot(vmin=vmin,#np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
               vmax=vmax,#np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)),
               cmap=cmap)
    plt.title('Coronal hole $EM_{'+wlen+'}$')
    plt.colorbar()
    plt.savefig(join(plotshome, 'ch_{:%Y-%m-%dT%H%M}{}_emthree').format(parse(chmap.date), thing, wlen))
    plt.close()
    print 'EM_three', label, np.nanmean(emmap.data)

    fig = plt.figure(figsize=(32, 24))
    emmap = Map(alldata/6, chmap.meta)
    emmap.plot(vmin=vmin,#np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
               vmax=vmax,#np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)),
               cmap=cmap)
    plt.title('Coronal hole $EM_{'+wlen+'}$')
    plt.colorbar()
    plt.savefig(join(plotshome, 'ch_{:%Y-%m-%dT%H%M}{}_emall').format(parse(chmap.date), thing, wlen))
    plt.close()

    plt.sca(hax)
    hist, bins = np.histogram(chmap.data.flatten(), range=(5.6, 7.0), bins=71)
    hist = (hist / float(chmap.size)) * 100
    plt.step(bins[:-1], hist, label=label, color=c)
    plt.xlabel('log(T)')
    plt.ylabel('% of image')

plt.legend()
plt.savefig(join(plotshome, 'ch-histograms'))
plt.close()
