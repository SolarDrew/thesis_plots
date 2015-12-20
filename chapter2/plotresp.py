from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
import numpy as np

tresp = readsav('/home/sm1ajl/CoronaTemps/aia_tresp')

tresp.resp94[:np.where(tresp.logt == 6.3)[0]] *= 6.7

fig = plt.figure(figsize=(16, 12))
plt.plot(tresp.logt, tresp.resp94, label='9.4nm')#, color='green')
plt.plot(tresp.logt, tresp.resp131, label='13.1nm')#, color='teal')
plt.plot(tresp.logt, tresp.resp171, label='17.1nm')#, color='gold')
plt.plot(tresp.logt, tresp.resp193, label='19.3nm')#, color='brown')
plt.plot(tresp.logt, tresp.resp211, label='21.1nm')#, color='purple')
plt.plot(tresp.logt, tresp.resp335, label='33.5nm')#, color='blue')
#plt.axvline(5.6, linestyle='--')
#plt.axvline(7.0, linestyle='--')
for r in [tresp.resp94, tresp.resp131, tresp.resp171, tresp.resp193, tresp.resp211, tresp.resp335]:
    print r.min(), r.max()
plt.yscale('log')
plt.xlim(4.5, 7.5)
plt.ylim(10**-28, 2*10**-24)
plt.legend(fontsize=20)
plt.title('AIA temperature response', fontsize=28)
plt.xlabel('log(T)', fontsize=24)
plt.ylabel('Temperature response [DN cm$^{5}$ s$^{-1}$ pix$^{-1}$]', fontsize=24)
plt.savefig('t_resps')
plt.close()
