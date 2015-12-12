from matplotlib import use, rc
use('agg')
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import sys
from os.path import expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as tm
import sunpy.map

