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
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/EISF_C6H12f', (u'netcdf',))
parameters['projection'] = None
parameters['q_vectors'] = 'SphLatt_2_50_2'
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/C6H12f_traject.nc'
parameters['weights'] = u'b_incoherent'

################################################################
# Setup and run the analysis                                   #
################################################################

eisf = REGISTRY['job']['eisf']()
eisf.run(parameters,status=True)