## Union components
The Union components are a bit different than other McStas components because they need to work together in order to function.
In this exercise you will add a sample environment to a simple instrument file called open_exercise.instr. The file already has a source, sample made with Union components and two detectors. The also file contains definitions of many materials, for example Aluminum and a Na2Ca3Al2F14 powder.

### TASK 1
* Add a Al sample environment to the instrument file using additional Union geometry components.
* Try to add small parts to the sample environment of other materials, it does not have to be realistic.
* Run simulations for a beam with energy of 3 meV, 4 meV and 5 meV.

### HINTS
* The higher priority is simulated when geometries overlap, the outer parts of a sample environment should have lower priority than the sample
* Use the 3D view to ensure that the geometry you build is correct
* You can make absorption images of your sample environment by making the source larger
* mcdoc Union
* mcdoc process
* Use logarithmic axes in plot by pressing L
* A solution is provided called open_exercise_solution.instr

### INTERPRETATION
Without any additions, the instrument returns a powder diffraction pattern. 
With an Al sample environment, the 3 different beam energies result in very different Al scattering. 
* With 3 meV no Al Bragg scattering is possible, no multiple scattering in cryostat
* With 4 meV only one Bragg reflection near backscattering is allowed
* With 5 meV two Bragg reflections are allowed


### TASK 2
Add Union logger components before the Union_master to record what happens in your sample environment. Run simulation for 3 meV, 4 meV and 5 meV beam energies to understand your previous results better.

### HINTS
* mcdoc Union_logger
* Loggers are placed in space, usual at the sample position
* The target_geometry option can select to only log in certain geometry, e.g. "sample_holder"
* The target_process option can select to only log a certain scattering process, e.g. "Al_powder"
* Need to use logarithmic axis in plotting, press "L"
* May need to use up to 1E8 rays to see necessary statistics

### INTERPRETATION
This interpretation section explains each logger result in the open_exercise_logger_solution.instr instrument. 
* scattering_time: Scattering as a function of time. No scattering happens before the beam reaches the sample environment, see the 5 early peaks that corresponds to entering two layers of the cryostat, the sample, and leaving the two layers. Then there is the multiple scattering decay that is highly dependent on the energy selected.
* scattering_time_sample: Scattering in sample as a function of time. The first main peak is the beam reaching the sample. At 5 meV energy two later peaks are observed, these correspond to the backscattered beam from when the beam tries to leave the cryostat.
* scattering_time_environment: Scattering in environment as a function of time. At 5 meV small peaks are seen when the beam bounces back and forth from backscattering in front and back.
* scattering_zx: Spatial distribution of scattering in zx plane. Beam entry, sample and beam exit are most intense. See powder scattering from beam entry and exit, but also from the powder sample. In general the scattered signals from the sample are thinner because the sample is smaller than the beam.
* scattering_xy: Spatial distribution of scattering in xy plane. The powder rings are clerkly visible from this view, and the sample holder can also be observed under the sample.
* scattering_zx: Spatial distribution of scattering in zy plane. Here the backscattering from entry and exit is more apparent, and can be distinguished from the scattering from the sample.
* reciprocal_space_zx: Scattering vector distribution in zx plane. The smaller circle on the right side correspond to the reach in reciprocal space of the initial beam, yet multiple scattering allows reaching more of reciprocal space. Each vertical line corresponds to a powder Bragg peak.



