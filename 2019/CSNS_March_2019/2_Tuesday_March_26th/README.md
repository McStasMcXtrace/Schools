This folder contains instrument files for use as solutions to exercises at the CSNS in March 2019.

The schedule of the day switches between talks about new components and practical sessions where
 these components are used in very simple instrument files.

## [A) Slits and collimators](1_slits_and_collimators)
	A1) Source, Slit, Slit, Monitor
	A2) Source, Collimator, Monitor

## [B) Guides](2_guides)
	B1) Source, Slit, Guide, Slit, Monitor (includes guide length as parameter)
	B2) Source, Guide, Monitor (eqvivalent but simpler)

## [C) Monochromator](4_monochromators)
	C1) Source, Guide, Monochromator, Monitor (includes 2theta calculation in initialize)

## [D) Choppers](5_choppers)
	D1) Source, Guide, Chopper, Guide, Monochromator, Monitor

## [E) Velocity Selector](6_velocity_selector)

These notes are for exercises held on Tuesday, but the students have
already had a practical session on sources and monitors on
monday. Assume sources and monitors are understood. Here we will use
```source_simple``` and ```PSD_monitor, Divergence_monitor```, ```L_monitor``` and ```DivPos_monitor```.


## A1) Slits
### TASK
Make a slit collimation system consisting of two slits.

### HINTS
* Slit component is called "Slit", for both circular and rectangular
* mcdoc ```Slit```
* Distances are important

### INTERPRETATION
Divergence diagram is not trivial, ray with max horizontal divergence can happen from all heights covered by the first slit, but with some vertical divergence a smaller area of the source can get through. Vice versa for vertical divergence. This makes the square divergence patter rotated 45 deg from the square slits. The acceptance diagram also shows the same story, and can be recreated by drawings.

## A2) Collimator
### TASK
Use a collimator instead of the slit system

### HINTS
* Simple collimator is called "Collimator_linear"
* ```mcdoc Collimator```
* Divergence often given in arc min, 60 arc min = 1 deg.

### INTERPRETATION
Somewhat like many small slit systems of larger vertical size than horizontal. Same argument for sloped divergence

## B) Guide
### TASK
Make a one segment guide and use a parameter to control its length

### HINTS
* Use Guide or Guide_gravity component
* ```mcdoc guide```
* Add parameter in DEFINE INSTRUMENT line between parenthesis
* m value important

### INTERPRETATION
Simple guide limits the divergence and intensity drop over distance

### EXTRA
Perform a scan over guide length to see intensity fall off

## C) Monochromator
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
The monochromator produces a low divergence beam. Contains correlation between divergence, position and energy due to Braggs law. Higher orders have different intensity due to guide transport for different wavelengths.

## D) Chopper
### TASK
Imagine if the higher orders were a problem, introduce a chopper into the guide to eliminate them

### HINTS
* The guide has to be split into two segments to accommodate the chopper
* Use DiskChopper, but set is_first to 0
* Ensure the distance between source and chopper is known
* Use declare and initialize to calculate appropriate delay
* McStas has a constant called K2V that converts from wavevector to velocity

### INTERPRETATION
Higher orders of scattering from monochromator should be removed.

## E) Velocity Selector
### TASK
Replace the higher order suppressing chopper with a velocity selector which only
allows the deisred wavelength through.

### HINTS
* ```mcdoc V_selector```
* You may (mostly) use the default values coming from D11@ILL
* You will have to make more room for the velocity selector than for he chopper
* Declare extra variable(s) to be able to dynamically set the wavelength
