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
parameters['frames'] = (0, 1000, 1)
parameters['grouping_level'] = u'atom'
parameters['instrument_resolution'] = ('gaussian', {'mu': 0.0, 'sigma': 0.1})
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/DISF_waterG', (u'netcdf',))
parameters['projection'] = None
parameters['q_vectors'] = 'SphLatt_3_21_2'
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/waterG_traject.nc'
parameters['weights'] = u'b_incoherent2'

################################################################
# Setup and run the analysis                                   #
################################################################

disf = REGISTRY['job']['disf']()
disf.run(parameters,status=True)