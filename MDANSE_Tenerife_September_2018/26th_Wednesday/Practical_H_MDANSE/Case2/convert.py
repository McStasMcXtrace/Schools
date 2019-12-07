#!/usr/bin/python

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE import REGISTRY

################################################################
# Job parameters                                               #
################################################################

parameters = {}
parameters['config_file'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/waterTIP4P2005_360mol.data'
parameters['mass_tolerance'] = 1e-05
parameters['n_steps'] = 2000
parameters['output_files'] = (u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/waterL_traject', (u'netcdf',))
parameters['time_step'] = 1.0
parameters['trajectory_file'] = u'/home/tofhr/gonzalezm/MDANSE2018/MDANSE_Case2/water.trj'

################################################################
# Setup and run the analysis                                   #
################################################################

lammps = REGISTRY['job']['lammps']()
lammps.run(parameters,status=True)
