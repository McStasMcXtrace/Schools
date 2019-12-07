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
parameters['instrument_resolution'] = ('gaussian', {'mu': 0.0, 'sigma': 0.1})
parameters['interpolation_order'] = u'1st order'
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/DOS_waterL', (u'netcdf',))
parameters['projection'] = None
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/waterL_traject.nc'
parameters['weights'] = u'equal'

################################################################
# Setup and run the analysis                                   #
################################################################

dos = REGISTRY['job']['dos']()
dos.run(parameters,status=True)