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
parameters['grouping_level'] = u'atom'
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/MSD_waterL', (u'netcdf',))
parameters['projection'] = None
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/waterL_traject.nc'
parameters['weights'] = u'equal'

################################################################
# Setup and run the analysis                                   #
################################################################

msd = REGISTRY['job']['msd']()
msd.run(parameters,status=True)