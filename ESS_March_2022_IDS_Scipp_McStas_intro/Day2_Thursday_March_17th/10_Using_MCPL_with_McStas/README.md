## MCPL practical

### Task a
* Use the SANS_spheres instrument from the example directory
* Insert an MCPL_output after last slit, geometrically at the slit
* Comment out the remainder of the instrument
* Run a series of runs with this new primary instrument at different statistics
* Create a "backend" instrument that uses MCPL_input as a source

### Task b
* Use ```ESS_butterfly_Guide_curved_test.instr```
* Decide for an ESS beamport of choice from https://public.esss.dk/users/willend/MCPL
* Use the information in
  https://confluence.esss.lu.se/display/MCSTAS/Using+MCPL+as+source+term+in+McStas
  to run comparative simulations of MCPL-based and analytical ESS
  source model

### Task c
* Investigate how MCPL is used in the model of ESS_BEER, combining
  primary-spectrometer from SIMRES and backend from McStas
