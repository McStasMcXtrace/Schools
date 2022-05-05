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
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/GR_waterG', (u'netcdf',))
parameters['r_values'] = (0.0, 1.0, 0.005)
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/waterG_traject.nc'
parameters['weights'] = u'equal'

################################################################
# Setup and run the analysis                                   #
################################################################

pdf = REGISTRY['job']['pdf']()
pdf.run(parameters,status=True)