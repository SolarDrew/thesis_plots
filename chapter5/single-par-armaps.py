from matplotlib import use, rc, patches, cm, _cm
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='25.0')
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import sys
import os
from os.path import join, expanduser, exists
CThome = join(expanduser('~'), 'CoronaTemps')
sys.path.append(CThome)
from temperature import TemperatureMap as tm
from sunpy.map import Map
from sunpy.time import parse_time as parse
from skimage import measure
from sunpy.net import hek

fs = 20
datahome = join('/fastdata', 'sm1ajl', 'thesis', 'data')
mapshome = datahome.replace('/data', '/maps')
plotshome = join(datahome.replace('/data', '/plots'), 'chapter5', 'ar')
if not exists(plotshome):
    os.makedirs(plotshome)

#dates = ['2011-01-22', '2011-01-20', '2011-02-01', '2011-02-19']
#coords = [([-200, 200], [300, 700]), ([-550, -250], [350, 650]),
#          ([-600, -200], [-450, -50]), ([-200, 200], [100, 500])]
dates = ['2011-01-22', '2011-02-01', '2011-02-19']
coords = [([-200, 200], [300, 700]), ([-600, -200], [-450, -50]), ([-200, 200], [100, 500])]
armaps = []

for date, coord in zip(dates, coords):
    d = parse(date)
    fdmap = Map(join(datahome, '{0:%Y/%m/%d}/171/*_{0:%Y?%m?%d}?{0:%H?%M?%S}*fits'.format(d)))
    if isinstance(fdmap, list): fdmap = fdmap[0]
    fig, ax = plt.subplots(figsize=(16, 12))
    fdmap.plot()
    corner = coord[0][0], coord[1][0]
    x, y = coord[0][1] - coord[0][0], coord[1][1] - coord[1][0]
    rect = patches.Rectangle(corner, x, y, color='white', fill=False)
    ax.add_artist(rect)
    plt.savefig(join(plotshome, 'fulldisk_ar_{:%Y-%m-%dT%H%M}').format(parse(date)))
    plt.close()
    thismap = tm(date, verbose=True,
                 data_dir=join(datahome, '{:%Y/%m/%d}'.format(d)),
                 maps_dir=join(mapshome, '{:%Y/%m/%d}'.format(d)))
    thismap.save()
    thismap = thismap.submap(*coord)
    armaps.append(thismap)

client = hek.HEKClient()
ar47 = client.query(hek.attrs.Time('2011-01-22', '2011-01-22'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11147')[2]
ar49 = client.query(hek.attrs.Time('2011-01-22', '2011-01-22'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11149')#[1]
ar50 = client.query(hek.attrs.Time('2011-02-01', '2011-02-01'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11150')[2]
ar61 = client.query(hek.attrs.Time('2011-02-19', '2011-02-19'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11161')#[1]
ar62 = client.query(hek.attrs.Time('2011-02-19', '2011-02-19'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11162')[0]

hfig, hax = plt.subplots(figsize=(16, 12))
for n, armap, arnum, stuff, label, c in [(0, armaps[0], 'AR11147 and AR11149', [ar47, ar49], 0, 0),
                                         (1, armaps[0], 'AR11147', [ar47], '$AR_1$', 'green'),
                                         (2, armaps[1], 'AR11150', [ar50], '$AR_2$', 'red'),
                                         (3, armaps[2], 'AR11161 and AR11162', [ar61, ar62], '$AR_3$', 'blue')]:
    if n == 0: continue
    fig = plt.figure(figsize=(32, 24))
    ax1 = fig.add_subplot(111)
    armap.plot()#vmin=5.8, vmax=6.5)
    #labels = ax1.get_xlabel(), ax1.get_ylabel()
    #ax1.set_xlabel(labels[0], fontsize=fs)
    #ax1.set_ylabel(labels[1], fontsize=fs)
    #for thing in stuff:
    #    rect = patches.Rectangle([thing, -1150], 700, 700, color='black', fill=False)
    #    ax1.add_artist(patches.Polygon(thing['hpc_bbox']))
    if len(stuff) > 1:
        plural = 's'
    else:
        plural = ''
    plt.title('Active region{} {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        plural, arnum, np.nanmin(armap.data), np.nanmean(armap.data), np.nanmax(armap.data)))#,
    #    fontsize=24)
    plt.colorbar()
    
    plt.savefig(join(plotshome, 'ar_{:%Y-%m-%dT%H%M}').format(parse(armap.date)))
    plt.close()

    vmin = 24.0
    vmax = 29.5
    alldata = np.zeros(shape=armap.shape)
    threedata = np.zeros(shape=(armap.shape))
    for wlen in ['94', '131', '171', '193', '211', '335']:
        fig = plt.figure(figsize=(32, 24))
        emmap = armap.calculate_em(wlen)
        alldata += emmap.data
        if wlen in ['171', '211', '193']: threedata += emmap.data
        emmap.plot(vmin=vmin,#np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
                   vmax=vmax)#np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)))
        plt.title('Active region $EM_{'+wlen+'}$')
        plt.colorbar()
        plt.savefig(join(plotshome, 'ar_{:%Y-%m-%dT%H%M}_em{}').format(parse(armap.date), wlen))
        plt.close()
    cmap = emmap.cmap

    fig = plt.figure(figsize=(32, 24))
    emmap = Map(threedata/3, armap.meta)
    emmap.plot(vmin=vmin,#np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
               vmax=vmax,#np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)),
               cmap=cmap)
    plt.title('Active region $EM_{three}$')
    plt.colorbar()
    plt.savefig(join(plotshome, 'ar_{:%Y-%m-%dT%H%M}_emthree').format(parse(armap.date), wlen))
    plt.close()
    print 'EM_three', label, np.nanmean(emmap.data)

    fig = plt.figure(figsize=(32, 24))
    emmap = Map(alldata/6, armap.meta)
    emmap.plot(vmin=vmin,#np.nanmean(emmap.data)-(2*np.nanstd(emmap.data)),
               vmax=vmax,#np.nanmean(emmap.data)+(2*np.nanstd(emmap.data)),
               cmap=cmap)
    plt.title('Active region $EM_{all}$')
    plt.colorbar()
    plt.savefig(join(plotshome, 'ar_{:%Y-%m-%dT%H%M}_emall').format(parse(armap.date), wlen))
    plt.close()

    plt.sca(hax)
    hist, bins = np.histogram(armap.data.flatten(), range=(5.6, 7.0), bins=71)
    hist = (hist / float(armap.size)) * 100
    plt.step(bins[:-1], hist, label=label, color=c)
    plt.xlabel('log(T)')
    plt.ylabel('% of image')

plt.legend()
plt.savefig(join(plotshome, 'ar-histograms'))
plt.close()
