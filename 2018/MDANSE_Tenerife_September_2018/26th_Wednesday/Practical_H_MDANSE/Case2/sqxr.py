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
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/SQXR_waterL', (u'netcdf',))
parameters['q_values'] = (2.0, 150.0, 1.0)
parameters['r_values'] = (0.0, 1.0, 0.005)
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/waterL_traject.nc'

################################################################
# Setup and run the analysis                                   #
################################################################

xssf = REGISTRY['job']['xssf']()
xssf.run(parameters,status=True)