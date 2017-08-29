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
rc('font', size='25.0')
import matplotlib.pyplot as plt
from matplotlib import patches
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
from sys import argv, exit

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
temps = np.arange(4.6, 7.405, 0.01)
widths = np.array([0.01, 0.1, 0.5]) # Just copying Aschwanden's range here
#widths = np.array([0.1]) # Just copying Aschwanden's range here
heights = 10 ** np.arange(18, 37, 0.1)
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

# Load AIA response functions
resp = load_temp_responses()

# Load unnessecary map for its metadata
voidmap = Map(sunpy.AIA_171_IMAGE)
mapmeta = voidmap.meta

e_model = np.memmap(filename=path.join('/fastdata', 'sm1ajl', 'synth_emiss_{}pars').format(n_pars),
                    dtype='float32', mode='r', shape=(141, 6))

# Run synthetic data through 1param tempmap method
for w, wid in enumerate(widths):
    for wl, wlength  in enumerate(['94', '131', '171', '193', '211', '335']):
        emiss = Map(emission[wl, :, w, :].copy(), mapmeta)
        emiss.cmap = sunpy.cm.get_cmap('sdoaia{}'.format(wlength))
        emiss.meta['naxis1'] = emiss.shape[1]
        emiss.meta['naxis2'] = emiss.shape[0]
        emiss.meta['cdelt1'] = np.log10(heights[1]) - np.log10(heights[0])
        emiss.meta['cdelt2'] = temps[1] - temps[0]
        emiss.meta['crval1'] = np.log10(heights[0])
        emiss.meta['crval2'] = temps[0]
        emiss.meta['crpix1'] = 0.5
        emiss.meta['crpix2'] = 0.5
        if wlength == '94': wlength = '094'
        fits_dir = path.join(datahome, 'synthetic', wlength)
        emiss.save(path.join(fits_dir, 'model.fits'), clobber=True)

        plt.plot()
        """fig = plt.figure(figsize=(8, 12))
        emiss.plot(cmap='cubehelix', aspect='auto')
        plt.colorbar()
        plt.savefig(path.join(outdir, 'emission',
                    'emiss_wid={:.3f}_wlen={}'.format(wid, wlength).replace('.', '_')))"""
        emiss.data /= emission[2, :, w, :]

    images = [Map(emission[i, :, w, :], mapmeta) for i in range(6)]
    if n_pars == 3:
        cmdargs = "mpiexec -n 10 python {} model {} {} {} {} {} {} {}".format(
            tmap_script, n_pars, datahome, None, None, True, True, True).split()
    else:
        cmdargs = "python {} model {} {} {} {} {} {} {}".format(
            tmap_script, n_pars, datahome, None, None, True, True, True).split()
    status = subp.call(cmdargs)
    newmap = TemperatureMap(fname=path.join(CThome, 'temporary.fits'))
    subp.call(["rm", path.join(CThome, 'temporary.fits')])
    data, meta = newmap.data, newmap.meta
    fitsmap = Map(newmap.goodness_of_fit.copy(), newmap.meta.copy())
    fitsmap = Map(np.log10(fitsmap.data), fitsmap.meta)
    fitsmap.data[np.where(newmap.goodness_of_fit == 0.0)] = 100000
    fitsmap.data[np.where(fitsmap.data == 100000)] = fitsmap.min()
    print '\n\n======', wid, fitsmap.min(), fitsmap.max(), '======\n\n'
    newmap.data = data

    truetemp = np.array(list(temps)*n_heights).reshape((n_heights, n_temps)).T
    diff = Map((abs(truetemp - data) / truetemp) * 100, newmap.meta.copy())
    if n_pars == 3:
        wdata = Map(newmap.dem_width, newmap.meta.copy())
        truew = np.ones(shape=(n_temps, n_heights)) * wid
        diffw = Map((abs(truew - wdata.data) / truew) * 100, newmap.meta.copy())
    
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

    # Define input temperature and EM to check, calculate the appropriate indices and get the output values
    #t_in = 5.15
    #t_in = 6.45
    t_in = 7.1
    tindex = int((t_in - temps.min()) / (temps[1] - temps[0]))
    em_in = 10 ** 25.0
    emindex = np.where(np.float32(heights) == em_in)[0]#int((em_in - heights.min()) / (heights[1] - heights[0]))
    print '$$$$', tindex, emindex
    t_out_all = newmap.data[:, emindex]
    t_out = newmap.data[tindex, emindex]
    em_out = emdata.data[tindex, emindex]
    if n_pars == 3:
        w_out = wdata.data[tindex, emindex]
    
    gof = fitsmap.data[tindex, emindex][0]
    gof_out = fitsmap.data[:, emindex]
    print tindex, temps[tindex], t_out
    print emindex, heights[emindex], 10 ** em_out

    """# Plot input v output log(T) for given input EM for single-parameter method
    fig = plt.figure(figsize=(14, 12))
    plt.plot(temps, t_out_all, color='black', linewidth=1.3)
    plt.axvline(5.6, linestyle='-.', color='blue')
    plt.axvline(7.0, linestyle='-.', color='blue')
    plt.axvline(t_in, linestyle='dotted', color='red')
    plt.axhline(t_out, linestyle='dotted', color='red')
    plt.plot(temps, temps, linestyle='--', color='green')
    plt.plot(t_in, t_out, 'x', color='black')
    plt.xlabel('Input log(T)')
    plt.ylabel('Output log(T)')
    plt.savefig(path.join(outdir, '{}par-temp-in-v-out-_wid={:.3f}'.format(n_pars, wid).replace('.', '_')))
    plt.close()

    # Plot goodness-of-fit for all solutions at given input EM
    fig = plt.figure(figsize=(14, 12))
    plt.plot(temps, gof_out, color='black', linewidth=1.3)
    plt.axvline(5.6, linestyle='-.', color='blue')
    plt.axvline(7.0, linestyle='-.', color='blue')
    plt.xlabel('Input log(T)')
    plt.ylabel('Solution GoF')
    plt.savefig(path.join(outdir, '{}par-temp-in-v-gof-_wid={:.3f}'.format(n_pars, wid).replace('.', '_')))
    plt.close()

    # Plot input DEM and calculated output for comparison
    fig, ax = plt.subplots(figsize=(21, 18))
    print t_in, wid, np.log10(em_in)
    plt.plot(logt, gaussian(logt, t_in, wid, np.log10(em_in)), color='green', label='Input DEM')
    if n_pars == 1:
        plt.plot(logt, gaussian(logt, t_out, 0.1, em_out), color='red', linestyle='--', label='Solution DEM')
    else:
        plt.plot(logt, gaussian(logt, t_out, w_out, em_out), color='red', linestyle='--', label='Solution DEM')
    plt.legend()
    plt.xlim(3, 9)
    plt.xlabel('log(T)')
    plt.ylabel('log(EM)')
    plt.text(3.3, 23, 'log(g) = {:.3f}'.format(gof))
    plt.savefig(path.join(outdir, '{}par-tin={:.2f}_wid={:.3f}_wlen={}'.format(n_pars, t_in, wid, emwlen).replace('.', '_')))
    plt.close()"""

    lwr = np.where(temps >= 5.6)[0][0]-1
    upr = np.where(temps <= 7.0)[0][-1]+1
    thist = temps[lwr:upr]
    """fig, ax = plt.subplots(figsize=(21, 18))
    e171 = Map(path.join(datahome, 'synthetic', '171', 'model.fits')).data[tindex, emindex]
    totalerr = 0
    for wl, wlength  in enumerate(['094', '131', '171', '193', '211', '335']):
        emiss_in = Map(path.join(datahome, 'synthetic', wlength, 'model.fits'))
        e_in = (emiss.data[:, emindex])# / e171)#[lwr:upr, 0]
        print wl, wlength, e_in.mean()
        #plt.plot(temps, np.log10(e/e171), label=str(int(wlength)/10.0)+'nm')
        plt.plot(temps, np.log10(e_in), label=str(int(wlength)/10.0)+'nm')
        #plt.plot(thist, np.log10(abs(e_in - e_model[:, wl])), label=str(int(wlength)/10.0)+'nm',
        #         linestyle='--')
        #totalerr = totalerr + abs(e_in - e_model[:, wl])
    #plt.plot(thist, np.log10(totalerr / 6.0), color='black', linewidth=1.3, label='GoF')
    plt.axvline(t_in, color='green')
    plt.axvline(t_out, linestyle='--', color='red')
    plt.legend(loc='lower right')
    plt.savefig(path.join(outdir, 'emission/{}par-tin={:.2f}_wid={:.3f}'.format(n_pars, t_in, wid, emwlen).replace('.', '_')))
    plt.close()"""

    allfits = np.load('full-fit-values.npy')
    if n_pars == 1:
        print allfits.shape
        thisfit = allfits[tindex, emindex][0]
    else:
        method_parvals = np.load('parvals.npy')
        #for i in method_parvals: print 'iiii', i
        print '\n\n****0', allfits.shape, allfits.min(), allfits.max()
        i = allfits.argmin()
        print '****1', i, method_parvals[i][0], method_parvals[i][1], np.log10(method_parvals[i][2])
        iw = int((method_parvals[i][1] - 0.1) / 0.1)
        iEM = int((np.log10(method_parvals[i][2]) - 20) / 0.1)
        i0 = iEM + (iw*151)
        print '****2', i0#, i1
        thisfit = allfits[i0::7*151]
        print '****3', thisfit.shape, thisfit.min(), thisfit.max()

    fig, ax = plt.subplots(figsize=(21, 18))
    plt.plot(thist, np.log10(thisfit))
    #plt.plot(thist, thisfit)
    plt.axvline(t_in, color='green', label='Input temperature')
    plt.axvline(t_out, linestyle='--', color='red', label='Solution temperature')
    plt.xlabel('log(T)')
    plt.ylabel('Goodness-of-fit')
    plt.legend()
    plt.savefig(path.join(outdir, 'fit_{}par-tin={:.2f}_wid={:.3f}'.format(n_pars, t_in, wid, emwlen).replace('.', '_')))
    plt.close()

