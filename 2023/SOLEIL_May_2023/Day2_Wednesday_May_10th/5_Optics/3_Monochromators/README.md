# Exercise: Monochromator
One of the most important devices of a beamline is of course the monochromator.
In the exercise we will set up a monochromator to pick out a single (part of) an undulator peak.
We will be using an Si(111) monochromator as an example.

<img src="http://pd.chem.ucl.ac.uk/pdnn/inst2/singleb.gif">

## 1. Insert the monochromator.
McXtrace now includes several different (but similar) monochromator models. We will start out with a simple Bragg -crystal monochromator.
Start either with a completely new instrument file - or continue with (a copy of) the earlier instrument with the focusing optics.
Insert a `Bragg_crystal` (from optics), make sure that it is big enough to accept the beam. Otherwise you can accept the default values to insert a Si-crystal and use its \<hkl\>=111 reflection.

As we did for the mirrors in the previous exercise the crystal needs to be ROTATED to match the desired energy. Using Bragg's law *n&lambda;=2d sin(&theta;)* (and &lambda; = 12.39842 / E) we may find the angle we need. Some hints:
- It is useful to define the rotation angle &theta; as an input parameter. This will enable us to do rocking curve scans.
- The crystal unit cell side length is 5.4309 Å
- Make sure that the source you are using does sample the energy your crystal picks out.
- To improve the computation efficiency, you may restrict the energy range around the nominal value.

Again, as before, you need an `Arm` after the crystal to make the Z-axis lie on the beam axis of your instrument. Monitor the beam reflected by the crystal using e.g. an `E_monitor`. Do you get what you expect?

If you do not get any rays through to the detector of your instrument, one way to figure what is wrong is to start a visual tracing. This can give you some hints. Another useful tool is to put a `PSD_monitor_4PI` at the crystal location after the crystal. This monitor will tell you in what direction the rays are travelling.

Try to perform a rocking curve of the monochromator to verify that the peak-shape is as you expected.

## 3. Switch the crystal to become a curved one - simulating a heat-bump.
Introduce a radius of curvature of 20 m to the crystal.
This is a simple way of modelling a heat-bump on the crystal.

1. To measure the introduced extra divergence we shall of course need a measure of the divergence for a flat monochromator. Therefore insert a ```Divergence_monitor``` downstream of the monochromator. Make sure it is large enough to measure the full beam. Set either ```nv=1``` or ```nh=1``` to target 1D horizontal and vertical divergence monitors respectively. Run a simulation to get a baseline measure of the divergence. 
 
2. Replace the ```Bragg_crystal monochromator``` with a ```Bragg_crystal_bent``` instance. The model is built to accommodate curvature of the crystal hull as well as the underlying crystal planes. These do not necessarily have to be the same, but for simplicity we will in the following assume that they are. Insert sensible values for `y_b`, `z_c`, `lattice_y_b`, `lattice_z_c` for your monochromator, and verify that the divergence changes.
3. Perform a rocking curve also of the bent monochromator, and check that the peak shape changes.

4. Do a few simulations at various radii of curvature
:bulb: Hint: it may be clever to define the radius as an input parameter that can be scanned. Now - let's try to plot the divergence profiles in the same plot. The data recorded by McXtrace monitors is stored as flat ascii files, where header lines are marked as comments beginning with a '#'. Using your regular plotting tool (be it matlab, gnuplot or similar), plot the divergence in the same figure:  
If you did a scan plotting could be something like (with gnuplot):
```gnuplot
plot 'my_instrument_20191203_0123/0/divergence_monitor.dat' u 1:2 w lp
replot 'my_instrument_20191203_0123/1/divergence_monitor.dat' u 1:2 w lp
```
or with matlab or octave:
```matlab
b=textread('my_instrument_20191203_0123/0/divergence_monitor.dat','','commentstyle','shell');
c=textread('my_instrument_20191203_0123/1/divergence_monitor.dat','','commentstyle','shell');
plot(b(:,1),b(:,2),'b',c(:,1),c(:,2),'r');
```

## 4. Add another crystal to create a DCM
We will now add another crystal set on the same rotation arm to create a DCM-system.  The first crystal in our example diverts the beam "upwards". 

<img src="http://pd.chem.ucl.ac.uk/pdnn/inst2/doubleb.gif">

The simplest way to add another crystal downstream from the first is to:
1. Insert a second crystal at a distance RELATIVE to the extra Arm you put after the first crystal
2. ROTATE the second crystal by the same angle as the first one but negative around the X-axis.
3. Add another Arm AT the same position as the crystal to put the beam axis along the Z-axis.

We may note a few things:
- This puts each crystal on its own rotation axis.
- The crystals are automatically aligned such that the beam always hits the crystal centre.
- The second crystal is actually hit on its back from "below"

4. Try to target a more realistic system where the crystals are set on the same rotation stage (or indeed are channel cut as in `SOLEIL_ROCK`), by setting up a single rotation Arm and then positioning the crystals relative to that.
5. For a curved crystal hitting from "below" is not the same as from "above". Using Arms - try to make sure the second crystal is in fact hit from "above".

## Now switch the DCM to become a Laue-monochromator.
Change the Bragg-monochromator (reflection) into a Laue (transmission). This is a rather simple matter of replacing the crystal components with the contributed COMPONENT Laue_crystal_BC.
1. Insert at Laue_crystal_BC at the site of the first monochromator
2. The default orientation of the Laue crystal is the same as for the Bragg crystals. This means that to have the beam traverse the smallest dimension we have to ROTATE it an extra 90 deg, as ```ROTATE (-90 -MONO_ANGLE,0,0) RELATIVE PREVIOUS```
3. As always add an extra Arm to point the z-axis along the beam. :bulb: Hint: this needs to also "undo" the extra 90 deg. Verify that your monchromator works.
4. Add a second crystal to make a double Laue monochromator. Work about how to align this using ROTATED. As do the other models - the Laue model only considers a single reflection (the 111 in this case). This means we have to either realign the inner crystal coordinate system or ROTATE the entire crystal 180 deg. In the first case this may be done using the parameter ```alphaz=-1``` (the default is to have set to 1)
5. Verify that your beamline looks as you want it to look like with a visual trace.
