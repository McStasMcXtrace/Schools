import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# X-ray diffraction data (Soper)
data = np.loadtxt('FxQ_water_Soper.dat')
qexp = data[:,0] / 0.1
sq1 = data[:,1]
sq1_error = data[:,2]
sq2 = data[:,3]
sq2_error = data[:,4]
sq3 = data[:,5]
sq3_error = data[:,6]

plt.plot(qexp, sq1, color = 'g', linestyle = ':', label = 'XR, norm IV')
plt.plot(qexp, sq2, color = 'b', linestyle = ':', label = 'XR, norm III')
plt.plot(qexp, sq3, 'o', color = 'r', label = 'XR, norm II')
ax = plt.gca()
ax.errorbar(qexp, sq3, yerr=sq3_error, color='r')

# X-ray S(Q) computed by MDANSE
f = Dataset('SQXR_waterL.nc', mode='r')
x_values = f.variables['q'][:]
x_name   = f.variables['q'].name
x_units  = f.variables['q'].units
y_values = f.variables['xssf_total'][:]
y_name   = f.variables['xssf_total'].name
y_units  = f.variables['xssf_total'].units
f.close()

plt.plot(x_values, 0.5*(y_values-1), color = 'k', linewidth=2.0, label = 'MD (no Q dep on f(Q))')
plt.xlabel(x_name + ' / ' + x_units)
plt.ylabel(y_name + ' / ' + y_units)
plt.xlim((0,100))
plt.legend()
plt.show()
