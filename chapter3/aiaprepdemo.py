from matplotlib import use, rc
use('agg')
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.instr.aia import aiaprep

# Open image from AIA during an SDO cross maneuver.
aiamap = sunpy.map.Map('aiacalibim1.fits')

# Process the AIA image
prepmap = aiaprep(aiamap)

# Plot original image against the processed image
fig = plt.figure(figsize=(20, 15))
fig.add_subplot(1, 2, 1)
aiamap.plot()
aiamap.draw_grid()
plt.title('Original image')
fig.add_subplot(1, 2, 2)
prepmap.plot()
prepmap.draw_grid()
plt.title('Processed image')
plt.savefig('aiatest_skimafftr_1')
plt.close()
