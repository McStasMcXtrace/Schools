#!/usr/bin/python

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE import REGISTRY

################################################################
# Job parameters                                               #
################################################################

parameters = {}
parameters['atom_aliases'] = {'CA': 'C', 'HA': 'H'}
parameters['field_file'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/FIELD'
parameters['history_file'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/HISTORY'
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case5/C6H12f_traject', (u'netcdf',))
parameters['version'] = u'2'

################################################################
# Setup and run the analysis                                   #
################################################################

dl_poly = REGISTRY['job']['dl_poly']()
dl_poly.run(parameters,status=True)