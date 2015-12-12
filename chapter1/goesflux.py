import matplotlib
from matplotlib import use, rc
use('agg')
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import sunpy
from sunpy.time import TimeRange
from sunpy.lightcurve import GOESLightCurve as glc
import datetime

flux = glc.create(TimeRange('2011/02/13 00:00', '2011/02/14 00:00'))
flux2 = glc.create(TimeRange('2011/02/14 00:00', '2011/02/15 00:00'))
flux3 = glc.create(TimeRange('2011/02/15 00:00', '2011/02/16 00:00'))
flux.data = flux.data.append(flux2.data).append(flux3.data)

figure = plt.figure(figsize=(16, 12))
axes = plt.gca()

dates = matplotlib.dates.date2num(flux.data.index._mpl_repr())

#axes.plot_date(dates, flux.data['xrsa'], '-',
#               label='0.5--4.0 $\AA$', color='blue', lw=2)
axes.plot_date(dates, flux.data['xrsb'], '-',
               label='0.1-0.8 nm', color='red', lw=2)

axes.set_yscale("log")
axes.set_ylim(1e-9, 1e-2)
axes.set_title("GOES X-ray Flux")
axes.set_ylabel('Watts m$^{-2}$')
#axes.set_xlabel(datetime.datetime.isoformat(flux.data.index[0])[0:10])
axes.set_xlabel("Universal Time")

ax2 = axes.twinx()
ax2.set_yscale("log")
ax2.set_ylim(1e-9, 1e-2)
ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))

axes.yaxis.grid(True, 'major')
axes.xaxis.grid(True, 'major')
axes.legend()

# @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
formatter = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M')
axes.xaxis.set_major_formatter(formatter)

axes.fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M')
figure.autofmt_xdate()

plt.savefig('flare-flux')
plt.close()
