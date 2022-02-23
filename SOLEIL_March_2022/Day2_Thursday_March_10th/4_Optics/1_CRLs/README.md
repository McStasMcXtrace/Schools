# Exercise 1: A set of compound refractive lenses
In this exercise we will and build a very simple beamline
based on compound refractive lenses (CRLs) and examine the focusing properties
hereof. 

## 1. Insert a source.
1. Insert a source. Position your cursor in the TRACE section below the
Progress_bar component and *Insert -> Sources/Source_gaussian*. 
1. Fill in the desired parameters:  

sig_x   | sig_y | sigPr_x | sigPr_y | gauss | dist | E0  | dE 
  :---:  | :---: | --- | --- | --- | --- | --- | ---
 48.2e-6 m | 9.5e-6 m | 100e-6 rad | 4.3e-6 rad |        1 |     31.5 m | 23.32 keV |   1 keV

and insert the coordinates for the source model: ( 0,0,0 ), and which object it should
be placed in relation to: *Origin*. Once you’re finished press OK to insert the source.
Note that we could just as easily have written this directly in the file using the text editor, or indeed any text-editor of
choice (Emacs, gedit etc...). McXtrace merely interprets the instrument file.
  
3. Try out the online component documentation by selecting Help/Component
 Library index. This should open a browser window with a listing of all
 available McXtrace components. From this web page you may access the online
 documentation for components and also the source code for the components
 should you wish to do so. Try to find the Source_gaussian model and what
 the parameters mean.

## 2. Monitor the inserted source model
1. Insert a 2D-monitor, and an energy spectrum monitor at the CRL-stack
   position 31.5 m downstream from the source. In our example this translates
   to a position of (0,0,31.5) in McXtrace coordinates. Use the components
   `PSD_monitor` and `E_monitor` to accomplish this. Contrary to other McXtrace
   components, a Monitor may reside on top of another component, if the
   restore_xray parameter is nonzero.
   Be sure to make the monitor large enough to catch the full beam footprint.
   Also, please set the optional filename parameter in each of the monitors.
   This is where your monitor data will be put. If none is given, the data are
   simply discarded and only the integrated intensity signal of the monitor is
   reported on the console. Run a simulation to get the beam footprint and spectrum at the
   CRL-stack position.
   
   To actually run the simulation just press Start. One may name an input
   directory in the “Output” box, if not McXtrace will generate a subdirectory
   based on the name of the beamline/instrument file and the current date and
   time. Once the simulation has run, you should be returned to the main
   window. To take a look at your newly generated data please press the “Plot”
   button. If you happen to have multiple CPUs/cores and MPI installed on your
   machine you may take advantage of this using the “Clustering” drop-down box
   and choose “MPI”.
   Is the result what you expected?

2. Optimize your simulation using directional sampling as you have done before today. Add the parameters
    `focus_xw=0.001`, `focus_yh=0.001`, and `dist=31.5` to your source parameters.
    This will cause McXtrace to only sample that part of the photon phase space
    which lies on trajectories that passing through a `focus_xw`  ×  `focus_yh` m^2
     aperture `dist` m downstream from the source. The weights
    of the photons are reduced to avoid biasing - this way we get correct
    intensity measures inside the sampling window.

## 3. Insert a CRL-stack.
The first stack should be a set of 16 parabolic Be lenses to set the focus at
10 m downstream of the source, the second stack a set of 16 Al lenses. Use
`Lens_parab` to set a Be CRL-stack 31.5 m downstream from the source, with the
following parameters:

Material |radius of curvature (at tip) | Number of lenses | Diameter | Thickness (at tip) 
--- | --- | --- | --- | ---
Be   | 200 &mu;m |   16 | 1 mm   | 50 &mu;m

## 4. Verify that the focus is where you expect

 1. Insert a 2D monitor at a variable position downstream of the source to
    verify that the focus is where you think it is. Make the z-position of the
    monitor a variable defined as an input parameter in the DEFINE INSTRUMENT
    line of your instrument file. This way you an avoid recompilation of your
    simulation every time you change the distance. Furthermore it allows you to
    perform pointwise scans of the parameter. In the appropriate boxes in the
    run-window ![scans](images/scan_L.png?raw=true "") you may input a range as min,max and a number of
    scan points.  
    TRICK: To really make it clear, you might replace the *Source_gaussian* with a *Source_div* (A 
    simple source model with a defined divergence) and set the divergence to 0. Thus you get a
    completely collimated source. Using what you know from earlier - can you do something
    like that with an EXTEND-block on the source?
    
 2. Use the lens stack together with a slit (Slit) as a low resolution
    monochromator, exploiting the achromatic focusing of the CRL. How would you
    set the slit? How would you verify the monochromatizing behaviour?

## 5. Visually trace photons through the beamline.

Use mxdisplay to try and visually track photons through the beamline model. In
the GUI this is done by selecting “Trace” instead of simulate in the run
window. You will need to zoom in a lot to see the CRLs clearly. You should be
able to see something like the figure below
[Visual Trace](images/trace.png)

## 6. Insert an Al stack of lenses

 1. Insert a second stack of CRLs with parameters at the same position as
    before:

Material | radius of curvature (at tip)  |   Number of lenses |   Diameter    |  Thickness (at tip)
--- | --- | --- | --- | ---
 Al   |     200 &mu;m            |       10|      1 mm      |    20 &mu;m

 To replace the CRL you could simply enclose the existing Lens_parab  COMPONENT with c-style comment markers:
 
 <code>
  /*
   COMPONENT crl = Lens_parab(
  ...
  )
 */
</code>

 and the either copy-paste the COMPONENT or `Insert` a new one using the second set of parameters.

 Is there a difference in performance between the CRL-stacks?
