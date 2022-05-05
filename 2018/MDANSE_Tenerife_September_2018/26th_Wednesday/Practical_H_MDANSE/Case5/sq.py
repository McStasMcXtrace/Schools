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
parameters['atom_transmutation'] = [[u'C6H12f_sel_H', u'h2']]
parameters['frames'] = (0, 1000, 1)
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/SQ_C6H12f', (u'netcdf',))
parameters['q_values'] = (3.0, 150.0, 1.0)
parameters['r_values'] = (0.0, 1.0, 0.005)
parameters['running_mode'] = ('monoprocessor',)
parameters['trajectory'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/C6H12f_traject.nc'
parameters['weights'] = u'b_coherent'

################################################################
# Setup and run the analysis                                   #
################################################################

ssf = REGISTRY['job']['ssf']()
ssf.run(parameters,status=True)