from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from matplotlib import patches
#from matplotlib import cm
import numpy as np
from temperature import TemperatureMap as tm
from sunpy.map import Map
from sunpy.time import parse_time as parse
from skimage import measure
from sunpy.net import hek
from os.path import join, expanduser

fs = 20
CThome = join(expanduser('~'), 'CoronaTemps')

qsmap1 = tm('2011-01-28', verbose=True, n_params=3,
            data_dir=join(CThome, 'data'),
            maps_dir=CThome)
#qsmap1 = tm('2011-02-01'))
qsmap1.save()
qsmap2 = tm('2011-02-08', verbose=True, n_params=3,
            data_dir=join(CThome, 'data'),
            maps_dir=CThome)
qsmap2.save()
qsmap3 = tm('2011-02-21', verbose=True, n_params=3,
            data_dir=join(CThome, 'data'),
            maps_dir=CThome)
qsmap3.save()
#fig = plt.figure(figsize=(24, 18))
fig = plt.figure(figsize=(16, 12))
for qsmap in [qsmap1, qsmap2, qsmap3]:
    #qsmap.save()
    #qsmap.data[np.where(qsmap.data == 0.0)] = np.NaN"""
    """try:
        centre_x = qsmap.reference_pixel['x']
        centre_y = qsmap.reference_pixel['y']
        x_grid, y_grid = np.mgrid[-centre_x:centre_x-1, -centre_y:centre_y-1]
        r_grid = np.sqrt((x_grid ** 2.0) + (y_grid ** 2.0))
        qsmap.data[r_grid > centre_x * 1.15] = None
    except:
        pass
    qsmap.save()"""

    ax1 = fig.add_subplot(1, 1, 1)
    #vmin = np.nanmean(qsmap.data, dtype=np.float64) - (2.0*np.nanstd(qsmap.data, dtype=np.float64))
    #vmax = np.nanmean(qsmap.data, dtype=np.float64) + (2.0*np.nanstd(qsmap.data, dtype=np.float64))
    qsmap.plot()#vmin=vmin, vmax=vmax, cmap='coolwarm')
    #[-400, 300], [-1150, -450]
    labels = ax1.get_xlabel(), ax1.get_ylabel()
    ax1.set_xlabel(labels[0], fontsize=fs)
    ax1.set_ylabel(labels[1], fontsize=fs)
    rect = patches.Rectangle([-600, -450], 400, 400, color='black', fill=False)
    rect2 = patches.Rectangle([-400, -1150], 700, 700, color='black', fill=False)
    ax1.add_artist(rect)
    ax1.add_artist(rect2)
    plt.title('Full disk corona {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        qsmap.date, np.nanmin(qsmap.data), np.nanmean(qsmap.data, dtype='float64'), np.nanmax(qsmap.data)),
        fontsize=24)
    plt.colorbar()
    plt.savefig('fulldisk_{:%Y-%m-%dT%H%M}'.format(parse(qsmap.date)))
    plt.close()

qsmap1 = qsmap1.submap([200, 700], [-700, -200])
qsmap2 = qsmap2.submap([-400, 100], [-400, 100])
qsmap3 = qsmap3.submap([-500, 0], [-500, 0])
print (qsmap3.data[qsmap3.data == 5.6]).size
qsmap3.data[qsmap3.data == 5.6] = np.NaN

for qsmap in [qsmap1, qsmap2, qsmap3]:
    print np.nanmin(qsmap.data), np.nanmean(qsmap.data), np.nanmax(qsmap.data)
    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(1, 1, 1)
    qsmap.plot()#vmin=5.95, vmax=6.35, cmap='coolwarm')
    labels = ax1.get_xlabel(), ax1.get_ylabel()
    ax1.set_xlabel(labels[0], fontsize=fs)
    ax1.set_ylabel(labels[1], fontsize=fs)
    plt.title('Quiet sun region {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        qsmap.date, np.nanmin(qsmap.data), np.nanmean(qsmap.data), np.nanmax(qsmap.data)),
        fontsize=24)
    plt.colorbar()
    
    plt.savefig('qs_{:%Y-%m-%dT%H%M}'.format(parse(qsmap.date)))
    plt.close()

"""armap1 = tm('2011-01-22').submap([-200, 200], [300, 700])#submap([-550, -250], [350, 650])
armap2 = tm('2011-01-20').submap([-550, -250], [350, 650])
armap3 = tm('2011-02-01').submap([-600, -200], [-450, -50])
#armap3 = tm('2011-02-19').submap([600, 1000], [-550, -150])
armap4 = tm('2011-02-19').submap([-200, 200], [100, 500])
coords = [([-200, 200], [300, 700]), ([-550, -250], [350, 650]),
          ([-600, -200], [-450, -50]), ([-200, 200], [100, 500])]
client = hek.HEKClient()
ar47 = client.query(hek.attrs.Time('2011-01-22', '2011-01-22'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11147')[2]
ar49 = client.query(hek.attrs.Time('2011-01-22', '2011-01-22'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11149')[1]
ar50 = client.query(hek.attrs.Time('2011-02-01', '2011-02-01'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11150')[2]
ar61 = client.query(hek.attrs.Time('2011-02-19', '2011-02-19'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11161')[1]
ar62 = client.query(hek.attrs.Time('2011-02-19', '2011-02-19'),
                    hek.attrs.EventType('AR'),
                    hek.attrs.AR.NOAANum == '11162')[0]

for n, armap, arnum, stuff in [(0, armap1, 'AR11147 and AR11149', [ar47, ar49]),
                               (1, armap2, 'AR11147', [ar47]),
                               (2, armap3, 'AR11150', [ar50]),
                               (3, armap4, 'AR11161 and AR11162', [ar61, ar62])]:
    map171 = Map('/media/huw/SDO_data/171/{:%Y/%m/%d/}*t00_00_*fits'\
        .format(parse(armap.date)))
    if isinstance(map171, list): map171 = map171[0]
    map171.data /= map171.exposure_time
"""    """for wl in [94, 131, 171, 193, 211, 335]:
        wlmap = Map('/media/huw/SDO_data/{}/{:%Y/%m/%d/}*t00_00_*fits'\
            .format(wl, parse(armap.date)))
        if isinstance(wlmap, list):
            wlmap = wlmap[0]
        wlmap.data /= wlmap.exposure_time
        wlmap.data = wlmap.data / map171.data.copy()
        print np.nanmin(wlmap), np.nanmax(wlmap)
        fig = plt.figure(figsize=(32, 18))
        fig.add_subplot(1, 2, 1)
        armap.plot()
        fig.add_subplot(1, 2, 2)
        wlmap = wlmap.submap(*coords[n])
        wlmap.plot(vmin=np.nanmean(wlmap.data)-(np.nanstd(wlmap.data)*2),
                   vmax=np.nanmean(wlmap.data)+(np.nanstd(wlmap.data)*2))#vmin=0, vmax=1500)
        #plt.colorbar(orientation='horizontal')
        plt.savefig('{}_{}'.format(parse(armap.date), wl))
        plt.close()"""
"""    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(111)
    #scale_cen = armap.mean()
    #scale_wid = 2 * armap.std()
    #armap.plot(vmin=scale_cen-scale_wid, vmax=scale_cen+scale_wid, cmap='coolwarm')
    #armap.plot(vmin=np.nanmin(armap.data), vmax=np.nanmax(armap.data), cmap='coolwarm')
    armap.plot(vmin=6.0, vmax=6.6, cmap='coolwarm')
    labels = ax1.get_xlabel(), ax1.get_ylabel()
    ax1.set_xlabel(labels[0], fontsize=fs)
    ax1.set_ylabel(labels[1], fontsize=fs)
    #for thing in stuff:
    #    rect = patches.Rectangle([thing, -1150], 700, 700, color='black', fill=False)
    #    ax1.add_artist(patches.Polygon(thing['hpc_bbox']))
    if armap in [armap1, armap4]:
        plural = 's'
    else:
        plural = ''
    plt.title('Active region{} {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        plural, arnum, np.nanmin(armap.data), np.nanmean(armap.data), np.nanmax(armap.data)),
        fontsize=24)
    plt.colorbar()
    
    plt.savefig('ar_{:%Y-%m-%dT%H%M}'.format(parse(armap.date)))
    plt.close()"""


qsmap = tm('2011-02-01')
lr = qsmap.superpixel((16, 16), method='average')
qsmap1 = qsmap.submap([-500, 500], [0, 1000])
qsmap2 = qsmap.submap([-400, 300], [-1150, -450])
lr1 = lr.submap([-500, 500], [0, 1000])
lr2 = lr.submap([-400, 300], [-1150, -450])
qsmap3 = tm('2011-02-14')
lr3 = qsmap3.superpixel((16, 16), method='average').submap([-400, 400], [-1200, -400])
qsmap3 = qsmap3.submap([-400, 400], [-1200, -400])
#qsline2 = qsmap2.submap([0, 0.6], [-1150, -1150+(700*0.375)])
qsline3 = qsmap3.submap([0, 0.6], [-1200, -1200+(800*0.375)])
#y2 = qsline2.pixel_to_data()[0][0] - 1150
#y2 = y2[1::3]
y3 = qsline3.pixel_to_data()[0][0] - 1200
y3 = y3[1::3]
#temps2 = [np.mean(qsline2.data[i-1:i+2, 0]) 
#                for i in range(1, len(qsline2.data), 3)]
temps3 = [np.mean(qsline3.data[i-1:i+2, 0])
                for i in range(1, len(qsline3.data), 3)]
#temps = qsline.data[:, 0]
fig = plt.figure(figsize=(16, 12))
#fig.add_subplot(2, 1, 1)
plt.title('Coronal temperatures above coronal hole', fontsize=24)
"""plt.plot(y2, temps2, color='black')#, label='Coronal hole 2')
plt.xlabel('X-position [arcsec]')
plt.axvline(-qsmap2.rsun_arcseconds, linestyle='--')
plt.xlim(-1200, -960)
fig.add_subplot(2, 1, 2)"""
plt.plot(y3, temps3, color='black')#'red', label='Coronal hole 3')
plt.axvline(-qsmap3.rsun_arcseconds, linestyle='--')
#plt.ylim(5.6, 6.15)
plt.xlim(-1200, -900)
plt.ylabel('log(T)')
plt.savefig('limb_temps_20110214')

#for qsmap, thing in [(qsmap1, 'a'), (qsmap2, 'b'), (qsmap3, '')]:
for qsmap, thing, lowres in [(qsmap1, 'a', lr1), (qsmap2, 'b', lr2), 
                             (qsmap3, '', lr3)]:
    #qsmap.data[qsmap.data < 5.8] = np.NaN
    scale_cen = qsmap.mean()
    scale_wid = 2 * qsmap.std()
    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(1, 1, 1)
    qsmap.plot(vmin=scale_cen-scale_wid, vmax=scale_cen+scale_wid, cmap='coolwarm')
    labels = ax1.get_xlabel(), ax1.get_ylabel()
    ax1.set_xlabel(labels[0], fontsize=fs)
    ax1.set_ylabel(labels[1], fontsize=fs)
    qsmap.draw_limb(color='black')
    #plt.axvline(x=0, ymin=0, ymax=0.375, color='blue')
    plt.title('Coronal hole {}\nMin: {:.2f}, mean: {:.2f}, max: {:.2f}'.format(
        qsmap.date, np.nanmin(qsmap.data), np.nanmean(qsmap.data), np.nanmax(qsmap.data)),
        fontsize=24)
    plt.colorbar()
    
    """data = qsmap.data.flatten()
    fig.add_subplot(1, 2, 2)
    plt.hist(data.flatten(), bins=141, range=(5.6, 7.0), normed=True)
    plt.title('Histogram of temperatures')
    plt.xlabel('log(T)')
    plt.ylabel('% of image')"""
    #plt.savefig('ch_{:%Y-%m-%dT%H%M}{}'.format(parse(qsmap.date), thing))

    #lowres = tm('2011-02-14').superpixel((16, 16), method='average').submap([-400, 400], [-1200, -400])
    contours = measure.find_contours(lowres.data, 6.1)
    for contour in contours:
        contour *= lowres.scale['x']
        contour[:, 0] += lowres.yrange[0]
        contour[:, 1] += lowres.xrange[0]
        plt.plot(contour[:, 1], contour[:, 0], color='green')
        plt.xlim(*qsmap.xrange)
        plt.ylim(*qsmap.yrange)
    
    plt.savefig('ch_{:%Y-%m-%dT%H%M}{}_cropped'.format(parse(qsmap.date), thing))
    plt.close()
