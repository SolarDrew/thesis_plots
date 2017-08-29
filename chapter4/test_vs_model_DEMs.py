# -*- coding: utf-8 -*-
"""
Script to produce synthetic AIA data based on arbitrary model DEMs and test the
results of the tempmap code against the model.

Created on Mon Jul 28 16:34:28 2014

@author: Drew Leonard
"""

import numpy as np
from matplotlib import use, rc
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
rc('font', size='22.0')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import cm
import sunpy
from sunpy.map import Map
from os import path, makedirs
import sys
CThome = path.join(path.expanduser('~'), 'CoronaTemps')
sys.path.append(CThome)
from temperature import TemperatureMap
from utils import gaussian, load_temp_responses
import subprocess32 as subp
from itertools import product
from skimage import measure
from sys import argv

# Decide whether to assess single-parameter or full-Gaussian method
n_pars = int(argv[1])

# Define which wavelength to use for EM estimation with 1-parameter TMaps
emwlen = str(argv[2])

#outdir = path.join(CThome, 'validation', '{}pars'.format(n_pars))
outdir = path.join('/fastdata', 'sm1ajl', 'thesis', 'plots', 'chapter4')
datahome = path.join('/fastdata', 'sm1ajl', 'thesis', 'data')
tmap_script = path.join(CThome, 'create_tempmap.py')
if not path.exists(outdir): makedirs(outdir)

# Define parameter ranges
temps = np.arange(4.6, 7.405, 0.01)#0.005)
widths = np.array([0.01, 0.1, 0.5])#np.arange(0.01, 0.605, 0.005) # Just copying Aschwanden's range here
heights = 10 ** np.arange(18, 37, 0.1)#0.05)
n_temps = len(temps)
n_widths = len(widths)
n_heights = len(heights)
parvals = np.array([i for i in product(temps, widths, heights)])
n_vals = n_temps * n_widths * n_heights

# Create model DEMs and synthetic emission
emission = np.zeros((6, n_temps, n_widths, n_heights))
logt = np.arange(0, 15.05, 0.05)
resp = load_temp_responses()
delta_t = logt[1] - logt[0]
for p, params in enumerate(parvals):
    dem = gaussian(logt, *params)
    f = resp * dem
    t = np.where(temps == params[0])[0][0]
    w = np.where(widths == params[1])[0][0]
    h = np.where(heights == params[2])[0][0]
    emission[:, t, w, h] = np.sum(f, axis=1) * delta_t

#emission = emission / emission[2, :, :, :]
#print '----', emission[2, :, :, :].min(), emission[2, :, :, :].max()

# Load AIA response functions
resp = load_temp_responses()

# Load unnessecary map for its metadata
voidmap = Map(sunpy.AIA_171_IMAGE)
mapmeta = voidmap.meta

# Run synthetic data through 1param tempmap method
for w, wid in enumerate(widths):#heights):
    print '\nWidth:', wid
    #fig = plt.figure(figsize=(30, 12))
    for wl, wlength  in enumerate(['94', '131', '171', '193', '211', '335']):
        #emiss = Map(emission[wl, :, :, w], mapmeta)
        emiss = Map(emission[wl, :, w, :].copy(), mapmeta)
        emiss.cmap = sunpy.cm.get_cmap('sdoaia{}'.format(wlength))
        emiss.meta['naxis1'] = emiss.shape[1]
        emiss.meta['naxis2'] = emiss.shape[0]
        #emiss.meta['cdelt1'] = widths[1] - widths[0]
        emiss.meta['cdelt1'] = np.log10(heights[1]) - np.log10(heights[0])
        emiss.meta['cdelt2'] = temps[1] - temps[0]
        #emiss.meta['crval1'] = widths[0]
        emiss.meta['crval1'] = np.log10(heights[0])
        emiss.meta['crval2'] = temps[0]
        emiss.meta['crpix1'] = 0.5
        emiss.meta['crpix2'] = 0.5
        if wlength == '94': wlength = '094'
        fits_dir = path.join(datahome, 'synthetic', wlength)
        if not path.exists(fits_dir): print makedirs(fits_dir)
        emiss.save(path.join(fits_dir, 'model.fits'), clobber=True)
        emiss.data /= emission[2, :, w, :]
        """ax = fig.add_subplot(1, 6, wl+1)
        emiss.plot(aspect='auto', vmin=emiss.min(), vmax=emiss.max())
        plt.title('{}'.format(wlength))
        plt.xlabel('Input EM')
        plt.ylabel('Input log(T)')
        plt.colorbar()
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
        plt.axhline(5.6, color='white')
        plt.axhline(7.0, color='white')
    #plt.savefig(path.join(outdir, 'model_emission_h={}'.format(np.log10(wid)).replace('.', '_')))
    plt.savefig(path.join(outdir, 'model_emission_w={}'.format(wid).replace('.', '_')))
    plt.close()"""
    
    #images = [Map(emission[i, :, :, w], mapmeta) for i in range(6)]
    images = [Map(emission[i, :, w, :], mapmeta) for i in range(6)]
    if n_pars == 3:
        cmdargs = "mpiexec -n 10 python {} model {} {} {} {} {} {} {}".format(
            tmap_script, n_pars, datahome, None, None, True, True, False).split()
    else:
        cmdargs = "python {} model {} {} {} {} {} {} {}".format(
            tmap_script, n_pars, datahome, None, None, True, True, False).split()
    status = subp.call(cmdargs)
    newmap = TemperatureMap(fname=path.join(CThome, 'temporary.fits'))
    subp.call(["rm", path.join(CThome, 'temporary.fits')])
    data, meta = newmap.data, newmap.meta
    fitsmap = Map(newmap.goodness_of_fit.copy(), newmap.meta.copy())
    fitsmap = Map(np.log10(fitsmap.data), fitsmap.meta)
    fitsmap.data[np.where(newmap.goodness_of_fit == 0.0)] = 100000
    fitsmap.data[np.where(fitsmap.data == 100000)] = fitsmap.min()
    newmap.data = data
    print '-------------MINMAX:-------------'#, newmap.min(), newmap.max(), newmap.shape, 
    #print newmap.data[newmap.data == 0].shape, '----------\n'
    print 'GoF', fitsmap.min(), fitsmap.mean(), fitsmap.max()
    print 'T_out', newmap.min(), newmap.mean(), newmap.max()

    #truetemp = np.array(list(temps)*n_widths).reshape((n_widths, n_temps)).T
    truetemp = np.array(list(temps)*n_heights).reshape((n_heights, n_temps)).T
    #print truetemp.shape, data.shape
    diff = Map((abs(truetemp - data) / truetemp) * 100, newmap.meta.copy())
    print 'T_diff', diff.min(), diff.mean(), diff.max()
    if n_pars == 3:
        wdata = Map(newmap.dem_width, newmap.meta.copy())
        truew = np.ones(shape=(n_temps, n_heights)) * wid
        #print 'truew', truew.min(), truew.mean(), truew.max()
        diffw = Map((abs(truew - wdata.data) / truew) * 100, newmap.meta.copy())
        print 'w_out', wdata.min(), wdata.mean(), wdata.max()
        print 'w_diff', diffw.min(), diffw.mean(), diffw.max()
    
    #print wid, newmap.xrange, newmap.yrange, newmap.scale
    #print wid, diff.xrange, diff.yrange, diff.scale
    #fig = plt.figure(figsize=(8, 12))
    fig = plt.figure(figsize=(24, 12))
    fig.add_subplot(1, 3, 1)
    newmap.plot(cmap='coolwarm', vmin=5.6, vmax=7.0, aspect='auto')
    plt.colorbar()
    plt.title('Solution log(T)')
    plt.ylabel('Input log(T)')
    plt.xlabel('Input EM')#width', fontsize=24)
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    
    ax = fig.add_subplot(1, 3, 2)
    #print 'diff', diff.min(), diff.max()
    #print np.nanmin(diff.data), np.nanmax(diff.data)
    diff.plot(cmap='RdYlGn_r', aspect='auto', vmin=0.0, vmax=30.0)
    plt.colorbar()
    plt.title('Difference from input (%)')
    plt.xlabel('Input EM')
    plt.ylabel('')
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    
    ax = fig.add_subplot(1, 3, 3)
    #print 'fits', fitsmap.min(), fitsmap.max()
    #print np.nanmin(fitsmap.data), np.nanmax(fitsmap.data)
    print '\n\n========', wid, fitsmap.min(), fitsmap.max(), '======'
    print '========', wid, np.nanmin(fitsmap.data), np.nanmax(fitsmap.data), '======\n\n'
    fitsmap.plot(cmap='cubehelix', aspect='auto',#)# vmax=0.1)
                 vmin=max(np.nanmean(fitsmap.data)-(2.0*(np.nanstd(fitsmap.data))), np.nanmin(fitsmap.data)),
                 vmax=min(np.nanmean(fitsmap.data)+(2.0*(np.nanstd(fitsmap.data))), np.nanmax(fitsmap.data)))
                 #vmin=np.nanmean(fitsmap.data)-(2.0*(np.nanstd(fitsmap.data))),
                 #vmax=np.nanmean(fitsmap.data)+(2.0*(np.nanstd(fitsmap.data))))
    plt.colorbar()
    plt.title('log(Goodness-of-fit)')
    plt.xlabel('Input EM')
    plt.ylabel('')
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    #plt.savefig(path.join(outdir, 'tempsolutions_em={}'.format(np.log10(wid)).replace('.', '_')))
    plt.savefig(path.join(outdir, '{}par-tempsolutions_wid={:.3f}'.format(n_pars, wid).replace('.', '_')))
    plt.close()

    if n_pars == 3:
        emdata = Map(newmap.emission_measure, newmap.meta.copy())
    else:
        if emwlen == 'three':
            total = np.zeros(newmap.shape)
            for wl in ['171', '193', '211']:
                emdata = newmap.calculate_em(wl, model=True)
                total += emdata.data
            emdata.data = total/3.0
        elif emwlen == 'all':
            total = np.zeros(newmap.shape)
            for wl in ['94', '131', '171', '193', '211', '335']:
                emdata = newmap.calculate_em(wl, model=True)
                total += emdata.data
            emdata.data = total/6.0
        else:
            emdata = newmap.calculate_em(emwlen, model=True)
    trueem = np.array(list(np.log10(heights))*n_temps).reshape(n_temps, n_heights)
    diffem = Map((abs(trueem - emdata.data) / trueem) * 100, newmap.meta.copy())
    fig = plt.figure(figsize=(24, 12))
    #fig = plt.figure(figsize=(8, 12))
    ax = fig.add_subplot(1, 3, 1)
    print 'em_out', emdata.min(), emdata.mean(), emdata.max()
    print 'em_diff', diffem.min(), diffem.mean(), diffem.max()
    emdata.plot(cmap='coolwarm', aspect='auto',
                vmin=9.0, vmax=36.0)
    #            vmin=emdata.min(), vmax=emdata.max())
    #            vmin=np.log10(heights[0]), vmax=np.log10(heights[-1]))
    """contours = measure.find_contours(emdata.data, np.log10(heights[0]))
    for contour in contours:
        contour[:, 0] *= emdata.scale['y']
        contour[:, 1] *= emdata.scale['x']
        contour[:, 0] += emdata.yrange[0]
        contour[:, 1] += emdata.xrange[0]
        plt.plot(contour[:, 1], contour[:, 0], color='blue')
        plt.xlim(*emdata.xrange)
        plt.ylim(*emdata.yrange)
    contours = measure.find_contours(emdata.data, np.log10(heights[-1]))
    for contour in contours:
        contour[:, 0] *= emdata.scale['y']
        contour[:, 1] *= emdata.scale['x']
        contour[:, 0] += emdata.yrange[0]
        contour[:, 1] += emdata.xrange[0]
        plt.plot(contour[:, 1], contour[:, 0], color='black')
        plt.xlim(*emdata.xrange)
        plt.ylim(*emdata.yrange)"""
    plt.colorbar()
    x, y = np.meshgrid(np.arange(emdata.xrange[0], emdata.xrange[1], emdata.scale['x']),
                       np.arange(emdata.yrange[0], emdata.yrange[1], emdata.scale['y']))
    contours = plt.contour(x, y, emdata.data, levels=np.arange(10.0, 40.0, 5.0), colors='grey')
    plt.clabel(contours, inline=0, fmt='%1.1f', fontsize=20, colors='k')
    if n_pars == 3:
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.title('Solution EM')
    plt.ylabel('Input log(T)')
    plt.xlabel('Input EM')#width', fontsize=24)
    if n_pars != 3:
        plt.savefig(path.join(outdir, '{}par-emsolutions_wid={:.3f}_wlen={}_t-only'.format(n_pars, wid, emwlen).replace('.', '_')))
    #plt.close()
    
    ax = fig.add_subplot(1, 3, 2)
    #fig, ax = plt.subplots(figsize=(8, 12))
    #fig, ax = plt.subplots(figsize=(24, 12))
    diffem.plot(cmap='RdYlGn_r', aspect='auto',
                vmin=0, vmax=55)
    #            vmin=diffem.min(), vmax=diffem.max())
    if n_pars == 3:
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('Difference from input (%)')
    plt.xlabel('Input EM')
    plt.ylabel('')
    if n_pars != 3:
        plt.savefig(path.join(outdir, '{}par-emsolutions_wid={:.3f}_wlen={}_d-only'.format(n_pars, wid, emwlen).replace('.', '_')))
    #plt.close()
    
    ax = fig.add_subplot(1, 3, 3)
    #fig, ax = plt.subplots(figsize=(8, 12))
    fitsmap.plot(cmap='cubehelix', aspect='auto')
    if n_pars == 3:
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('log(Goodness-of-fit)')
    plt.xlabel('Input EM')
    plt.ylabel('')

    if n_pars == 3:
        plt.savefig(path.join(outdir, '{}par-emsolutions_wid={:.3f}'.format(n_pars, wid).replace('.', '_')))
    else:
        plt.savefig(path.join(outdir, '{}par-emsolutions_wid={:.3f}_wlen={}'.format(n_pars, wid, emwlen).replace('.', '_')))
    plt.close()

    if n_pars == 3:
        fig = plt.figure(figsize=(24, 12))
        ax = fig.add_subplot(1, 3, 1)
        wdata.plot(cmap='coolwarm', vmin=widths[0], vmax=widths[-1], aspect='auto')
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
        plt.axhline(5.6, color='white')
        plt.axhline(7.0, color='white')
        plt.colorbar()
        plt.title('Solution width')
        plt.ylabel('Input log(T)')
        plt.xlabel('Input EM')#width', fontsize=24)
    
        ax = fig.add_subplot(1, 3, 2)
        diffw.plot(cmap='RdYlGn_r', vmin=diffw.min(), vmax=diffw.max(), aspect='auto')
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
        plt.axhline(5.6, color='white')
        plt.axhline(7.0, color='white')
        plt.colorbar()
        plt.title('Difference from input (%)', fontsize=28)
        plt.ylabel('')
        plt.xlabel('Input EM')

        ax = fig.add_subplot(1, 3, 3)
        fitsmap.plot(cmap='cubehelix', aspect='auto')
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
        plt.axhline(5.6, color='white')
        plt.axhline(7.0, color='white')
        plt.colorbar()
        plt.title('log(Goodness-of-fit)')
        plt.ylabel('')
        plt.xlabel('Input EM')
    
        plt.savefig(path.join(outdir, '{}par-widsolutions_wid={:.3f}'.format(n_pars, wid).replace('.', '_')))
        plt.close()

"""w = np.where((widths > 0.097)*(widths < 0.103))
dataslice = data[:, w].reshape(len(temps))
diffslice = diff[:, w].reshape(len(temps))
fitslice = fits[:, w].reshape(len(temps))

fig = plt.figure(figsize=(16, 12))
plt.plot(temps, dataslice)
plt.title('Solution log(T) at width=0.1', fontsize=28)
plt.xlabel('Input log(T)', fontsize=24)
plt.ylabel('Solution log(T)', fontsize=24)
plt.savefig('/home/drew/Dropbox/euroscipy/dataslice')
plt.close()

fig = plt.figure(figsize=(16, 12))
plt.plot(temps, diffslice)
plt.title('Difference from input at width=0.1', fontsize=28)
plt.xlabel('Input log(T)', fontsize=24)
plt.ylabel('Difference (%)', fontsize=24)
plt.savefig('/home/drew/Dropbox/euroscipy/diffslice')
plt.close()

fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(1, 1, 1)
plt.plot(temps, np.log10(fitslice))
plt.title('Goodness-of-fit at width=0.1', fontsize=28)
plt.xlabel('Input log(T)', fontsize=24)
plt.ylabel('log(Goodness-of-fit)', fontsize=24)
plt.savefig('/home/drew/Dropbox/euroscipy/fitslice')
plt.close()"""
