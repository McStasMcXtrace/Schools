import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

pdf = Dataset('GR_waterL.nc', mode='r')
x_values = pdf.variables['r'][:]
x_name   = pdf.variables['r'].name
x_units  = pdf.variables['r'].units
y_values = pdf.variables['pdf_inter_HH'][:]
y_name   = pdf.variables['pdf_inter_HH'].name
y_units  = pdf.variables['pdf_inter_HH'].units
pdf.close()

# Neutron diffraction data (Soper)
data = np.loadtxt('soper13_water_structure_review_rdf.dat')
rnm = data[:,0] * 0.1
gr = data[:,5]
gr_error = data[:,6]
plt.plot(rnm, gr, 'o', color = 'r', label = 'Soper')
ax = plt.gca()
ax.errorbar(rnm, gr, yerr=gr_error, color='r')

plt.plot(x_values, y_values, color = 'k', linewidth=2.0, label = 'MD')
ax = plt.gca()
plt.xlabel(x_name + ' / ' + x_units)
plt.ylabel(y_name + ' / ' + y_units)
plt.xlim((0,1))
plt.legend()
plt.show()
