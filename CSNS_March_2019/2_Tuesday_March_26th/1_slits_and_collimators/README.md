## Slits
### TASK
Make a slit collimation system consisting of two slits, monitoring how
the beam/phase space changes. Use both ```PSD_monitor```, ```Divergence_monitor``` and ```DivPos_monitor```.

### HINTS
* Slit component is called "Slit", for both circular and rectangular
* ```mcdoc Slit```
* Distances are important

### INTERPRETATION
Divergence diagram is not trivial, ray with max horizontal divergence can happen from all heights covered by the first slit, but with some vertical divergence a smaller area of the source can get through. Vice versa for vertical divergence. This makes the square divergence patter rotated 45 deg from the square slits. The acceptance diagram also shows the same story, and can be recreated by drawings.

## Collimator
### TASK
Use a collimator instead of the slit system

### HINTS
* Simple collimator is called "Collimator_linear"
* ```mcdoc Collimator```
* Divergence often given in arc min, 60 arc min = 1 deg.

### INTERPRETATION
Somewhat like many small slit systems of larger vertical size than horizontal. Same argument for sloped divergence
