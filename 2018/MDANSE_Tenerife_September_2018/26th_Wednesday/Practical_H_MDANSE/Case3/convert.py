#!/usr/bin/python

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE import REGISTRY

################################################################
# Job parameters                                               #
################################################################

parameters = {}
parameters['fold'] = False
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/waterG_traject', (u'netcdf',))
parameters['pdb_file'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/TIP4P2005f-N500.pdb'
parameters['xtc_file'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case3/TIP4P2005f-N500.xtc'

################################################################
# Setup and run the analysis                                   #
################################################################

gromacs = REGISTRY['job']['gromacs']()
gromacs.run(parameters,status=True)
