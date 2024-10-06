# A powder diffractometer

Based on the output of the earlier exercises, we will now assemble a
powder diffractometer.

## TASKS
* Find back your solution to the [monochromator exercise](../../01_Monday_October_7th/04_Neutron_optics_exercises/Exercise_monochromator/) or take one one of those [provided](../../01_Monday_October_7th/04_Neutron_optics_exercises/Exercise_monochromator/solution)
* Using mcdoc and the PDF information, add `PowderN` to the sample position. Sample geometry should be a cylindar of `radius=0.005` and `yheight=0.07`. Sample definition can be one of your choice, but a good candidate is `reflections=Na2Ca3Al2F14.laz`
* Next, add a banana-shaped detector (cyllindrical cut) by using `Monitor_nD` of radius 1.2 m and height 30 cm, measuring a diffractogram using `options="banana theta bins=640 limits=[-170 -10]"` - positioned ` AT (0,0,0) RELATIVE ` your sample
* Run a simulation with `1e7` neutron rays, hopefully you have arrived at something that resembles this as a diffraction pattern:
![Initial diffractogram](pics/diffractogram1.png)