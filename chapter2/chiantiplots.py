from matplotlib import use, rc_params
use('agg')
import chianti
import chianti.core as ch
from chianti import filters
import numpy as np
import matplotlib.pyplot as plt

temperature = [1.0e6, 2.0e6]
density = 1.e9
wvl = np.arange(50.0, 450.0, 0.05)
emeasure = 1.e+25
s = ch.spectrum(temperature, density, wvl,em=emeasure,
                #filter=(filters.gaussian, 0.2), em=emeasure,
                doContinuum=0, minAbund=1.e-4)
wvl = [w /10.0 for w in wvl]

print s.Spectrum.keys()
print type(s.Spectrum['intensity'])

fig = plt.figure()
plt.plot(wvl, s.Spectrum['intensity'][0]*emeasure)
plt.ylabel('Intensity')
plt.xlabel('Wavelength (nm)')
plt.savefig('spectra1e6')
plt.close()

fig = plt.figure()
plt.plot(wvl, s.Spectrum['intensity'][1]*emeasure)
plt.ylabel('Intensity')
plt.xlabel('Wavelength (nm)')
plt.savefig('spectra2e6')
plt.close()

s.gofnt()
