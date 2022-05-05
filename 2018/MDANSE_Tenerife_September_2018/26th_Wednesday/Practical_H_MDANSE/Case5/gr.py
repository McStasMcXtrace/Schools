#!/usr/bin/python

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE import REGISTRY

################################################################
# Job parameters                                               #
################################################################

parameters = {}
parameters['atom_selection'] = None
parameters['atom_transmutation'] = None
parameters['frames'] = (0, 1000, 1)
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/GR_C6H12f', (u'netcdf',))
parameters['r_values'] = (0.0, 1.0, 0.005)
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/C6H12f_traject.nc'
parameters['weights'] = u'equal'

################################################################
# Setup and run the analysis                                   #
################################################################

pdf = REGISTRY['job']['pdf']()
pdf.run(parameters,status=True)