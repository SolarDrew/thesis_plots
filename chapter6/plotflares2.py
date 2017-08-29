# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 12:49:28 2014

@author: drew
"""

from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from sunpy import wcs
from sunpy.map import Map
from sunpy.net import hek, vso
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import numpy as np
import datetime as dt
from repeat_tempmaps import repeat
from astropy import units as u


class DownloadError(Exception):
    def __init__(self, msg):
        self.msg = msg

start = parse('2011-02-10')
end = parse('2011-02-16')
#end = parse('2011-01-04')

client = hek.HEKClient()
flares = client.query(hek.attrs.Time(start, end),
                      hek.attrs.EventType('FL'))

flares = [fl for fl in flares if (fl['ar_noaanum'] > 11137 and 
                                  fl['ar_noaanum'] < 11184)]

fluxes = []
temps = []
coords = []

for flare in flares:
  try:
    flaretime = parse(flare['event_starttime'])
    starttime = flaretime-dt.timedelta(hours=0.5)
    timerange = tr(starttime, flaretime)

    fits_dir = '/media/huw/SDO_data/171/{:%Y/%m/%d}/'.format(starttime)
    filename = fits_dir + 'aia*171*{0:%Y?%m?%d}?{0:%H?%M?%S}*lev1?fits'.format(
        flaretime)
    try:
        temp_im = Map(filename)
    except:
        print 'Failed to open map with filename ', filename
        raise
    # Download data if not enough found
    if temp_im == []:
        vsocl = vso.VSOClient()
        wlen = u.Quantity(value=171, unit='Angstrom')
        qr = vsocl.query(vso.attrs.Time(flaretime,
                                        flaretime + dt.timedelta(seconds=12)),
                          vso.attrs.Wave(wlen, wlen),
                          vso.attrs.Instrument('aia'),
                          vso.attrs.Provider('JSOC'))
        res = vsocl.get(qr, path=fits_dir+'{file}', site='NSO', 
                        methods=['URL_FILE_Rice']).wait()
        if res == []:
            res = vsocl.get(qr, path=fits_dir+'{file}', site='NSO').wait()
        temp_im = Map(res)
        if isinstance(temp_im, list):
            temp_im = temp_im[0]
    print flaretime
    region = client.query(hek.attrs.EventType('AR'),
                          hek.attrs.Time(flaretime-dt.timedelta(minutes=5), 
                                         flaretime))
    if flare['ar_noaanum'] == 11149:
        print 'HEK is fucked, using 11147'
        flare['ar_noaanum'] = 11147
    region = [r for r in region if r['ar_noaanum'] == flare['ar_noaanum']]
    if isinstance(region, list):
        try:
            region = region[0]
        except IndexError:
            'FAIL'
            continue
    coords = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'], 
                                b0_deg=temp_im.heliographic_latitude, 
                                l0_deg=temp_im.carrington_longitude)
    
    means, fluxes, maxes, times = repeat(starttime, flaretime, timeres=1.0/60.0, 
                                         coords=coords, ar=flare['ar_noaanum'], 
                                         plotminmax=True)
    
    fig = plt.figure(figsize=(18, 14))
    print len(times), len(maxes)
    plt.plot(times, maxes)
    plt.axhline(np.mean(maxes))
    plt.title('AR {} - {} flare'.format(flare['ar_noaanum'], 
                                        flare['fl_goescls']))
    plt.savefig('/home/drew/Dropbox/SIPWork/flareplots/AR{}_fl{}.png'.format(
        flare['ar_noaanum'], flare['fl_goescls']))
    plt.close()
  except:
      print 'FAIL!!!', flare['ar_noaanum'], flare['fl_goescls']
      raise
