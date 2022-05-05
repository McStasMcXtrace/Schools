## Phonon scattering
This exercise uses a triple axis instrument that investigates a Ni crystal in the HK,L plane. The instrument parameters are QHK and QL for selecting a point in reciprocal space using reciprocal lattice units. The initial and final energy are selected with EI and EF (in meV).

For example one can access the 222 Bragg peak with QHK=2 QL=2 EI=50 EF=50

Other allowed Bragg peaks include 002 and 004.

The sample is built using Union components and consists of incoherent scattering, single crystal scattering and simple phonon scattering. The processes are built on the Single_crystal and Phonon_simple components.

### TASK
Using the triple axis instrument and scanning try to find important features of the crystal
* Scan a Bragg peak in energy direction e.g. QHK=2 QL=2 EI=80 EF=60,100
* Scan a Bragg peak in reciprocal space, e.g. QHK=2 QL=1.5,2.5 and EI=EF
* Scan a phonon branch in energy direction, e.g. QHK=0 QL=3.5 EI=80 EF=50,70
* Scan both phonon branches in reciprocal space, e.g. QHK=0 QL=3,5 EI=80 EF=60
* See if you can find any instrument effects in the data

### HINTS
* It is important to select reasonable number of rays and scan points. 5E6 rays is sufficient, and less possible.
* Try the 3D visualization
* The instrument has 4 monitors. A monitor before the sample showing the intensity, and 3 PSD monitors after the sample with the same position. One records all events, one records only single order scattering and the last only multiple scattering.

### INTERPRETATION
A 2D scan can be seen in the longitudinal_scan.png file in this folder. 


