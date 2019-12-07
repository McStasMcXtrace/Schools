#!/usr/bin/python

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE import REGISTRY

################################################################
# Job parameters                                               #
################################################################

parameters = {}
parameters['atom_selection'] = [u'waterG_selOHH']
parameters['atom_transmutation'] = None
parameters['frames'] = (0, 1000, 10)
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/SQXR_waterG', (u'netcdf',))
parameters['q_values'] = (3.0, 150.0, 1.0)
parameters['r_values'] = (0.0, 1.0, 0.005)
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/waterG_traject.nc'

################################################################
# Setup and run the analysis                                   #
################################################################

xssf = REGISTRY['job']['xssf']()
xssf.run(parameters,status=True)