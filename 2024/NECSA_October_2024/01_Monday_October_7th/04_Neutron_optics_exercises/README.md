## TASK 1 - DiskChopper
Using a broadband source, let's try to limit the wavelength band coming through the optics.
You may start from Exercise_chopper.instr


## HINTS
* The guide has to be split into two segments to accommodate the chopper
* Use DiskChopper, but set is_first to 0
* Ensure the distance between source and chopper is known
* You may (but you don't have to) use declare and initialize to calculate appropriate delay
* McStas has a constant called K2V that converts from wavevector to velocity

You _may_ cheat by looking at the solution file

### TASK - Velocity_selector
Start from the instrument file [Exercise_vselector_start.instr](Exercise_vselector_start.instr)
Replace the higher suppressing chopper with a velocity selector which only
allows the deisred wavelength through.

### HINTS
* ```mcdoc V_selector```
* You may (mostly) use the default values coming from D11@ILL
* You will have to make more room for the velocity selector than for he chopper
* Declare extra variable(s) to be able to dynamically set the wavelength

