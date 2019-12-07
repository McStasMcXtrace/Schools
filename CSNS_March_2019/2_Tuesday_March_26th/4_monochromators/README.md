## Monochromators
### TASK
* Insert a monochromator after the guide and define wavelength input parameter
* Use declare and initialize to set the monochromator angle
* Observe second order

### HINTS
* Use Monochromator_curved component as it supports higher orders, for one crystal: NH=1,NV=1
* ```mcdoc Monochromator```
* Q of Pyrolytic Graphite = 3.355. Remember lambda < 4*PI/Q
* Variables declared in declare section with c syntax, RAD2DEG is the constant for radians to deg conversion
* Calculations prior to raytracing performed in initialize section
* Use Arm to find 2 theta direction

### INTERPRETATION
The monochromator produce a low divergence beam. Contains correlation between divergence, position and energy due to Braggs law. Higher orders have different intensity due to guide transport for different wavelengths. 
