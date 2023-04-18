# McXtrace training: components

Lecture describing components
1. What objects are around
2. How are they actually written? What do they contain/do?
3. Exercise: customize a monitor, adding a crude efficiency

## What are McXtrace components

The [McXtrace](https://www.mcxtrace.org) components are the building blocks of a beam-line description. In practice, they follow a similar syntax as beam-lines (DEFINE, TRACE, ...) but also include C code to specify how rays are transformed by each component. More specifically, components are the place where the physics takes place.

You can get an overview of all components from the Help menu and Docs button in the MxGUI window, or by typing the command `mxdoc` in a Terminal window. A browser window appears. Clicking on the *comp* item provides its description (header) and source code (and its location in the URL).

You can also access the list of all available components on-line [HERE](https://www.mcxtrace.org/download/components/3.1/).

## Overview of Source components

### Synchrotron Insertion Devices
- [Bending magnet](https://www.mcxtrace.org/download/components/3.1/sources/Bending_magnet.html)  see B.D. Patterson, Am. J. Phys. 79, 1046 (2011)
- [Undulator](https://www.mcxtrace.org/download/components/3.1/sources/Undulator.html) see K.J. Kim, AIP, conf. proc., 184, 1989
- [Wiggler](https://www.mcxtrace.org/download/components/3.1/sources/Wiggler.html) see B.D. Patterson, Am. J. Phys. 79, 1046 (2011)

### Lab/ideal stuff
- [Laboratory X-ray tube](https://www.mcxtrace.org/download/components/3.1/sources/Source_lab.html)  (e.g. rotating anode)
- Ideal, point and Gaussian

### Interfaces with other software

By its modularity and open source design, McXtrace can communicate with other X-ray simulation software (R: read, W: write).

- [Spectra](https://www.mcxtrace.org/download/components/3.1/sources/Source_spectra.html) (R) <http://spectrax.org/spectra/>
- [Simplex](https://www.mcxtrace.org/download/components/3.1/sources/Source_simplex.html) (R) <http://spectrax.org/simplex/index.html >
- [Genesis](https://www.mcxtrace.org/download/components/3.1/sources/Source_genesis13.html) (R) <http://genesis.web.psi.ch/>
- [Shadow](https://www.mcxtrace.org/download/components/3.1/misc/Shadow_input.html) (RW) <https://github.com/oasys-kit/shadow3>
- [MCPL](https://www.mcxtrace.org/download/components/3.1/misc/MCPL_input.html) (GEANT4, PHITS, MCNP,SRW) (RW) <https://mctools.github.io/mcpl/>
- SRW (R) <https://github.com/ochubar/SRW> Our converter generates an MCPL exchange file from SRW. You may also look at [OASYS](https://www.aps.anl.gov/Science/Scientific-Software/OASYS) to export interchange formats.

## Optics components

There is a dedicated session about this topic, but we here list a few components that are available to describe beam-lines.

- [Bragg crystal](https://www.mcxtrace.org/download/components/3.1/optics/Bragg_crystal.html) (monochromator, incl. [bent](https://www.mcxtrace.org/download/components/3.1/optics/Bragg_crystal_bent.html))
- [Capillary](https://www.mcxtrace.org/download/components/3.1/optics/Capillary.html)
- [Filter](https://www.mcxtrace.org/download/components/3.1/optics/Filter.html) (absorption and refraction)
- [Lenses](https://www.mcxtrace.org/download/components/3.1/optics/Lens_simple.html)
- Mirrors ([flat](https://www.mcxtrace.org/download/components/3.1/optics/Mirror.html), [curved](https://www.mcxtrace.org/download/components/3.1/optics/Mirror_curved.html), [multi-layers](https://www.mcxtrace.org/download/components/3.1/optics/Multilayer_elliptic.html), [twin KB](https://www.mcxtrace.org/download/components/3.1/optics/TwinKB_ML.html)  multi-layer)
- [Zone plate](https://www.mcxtrace.org/download/components/3.1/optics/ZonePlate.html)
- [Grating](https://www.mcxtrace.org/download/components/3.1/contrib/Reflective_grating.html) (lamellar, blazed)
- [Slit](https://www.mcxtrace.org/download/components/3.1/optics/Slit.html), [beam-stop](https://www.mcxtrace.org/download/components/3.1/optics/Beamstop.html), ...

## Samples

Samples are essential to build so-called virtual beam-lines, and reproduce data that look like real experiments. There are dedicated sessions on this topic.

To have an overview of existing sample, have a lokk at the [official samples](https://www.mcxtrace.org/download/components/3.1/#samples) and [contributed samples](https://www.mcxtrace.org/download/components/3.1/#contrib).

Here is a list of the main components you may use:

- [SaSView models](https://www.mcxtrace.org/download/components/3.1/samples/SasView_model.html): 60 models from SasView, PDB, Nanodiscs, Liposomes, I(q), ... [Saxs_sphere](https://www.mcxtrace.org/download/components/3.1/samples/Saxs_spheres.html), [SAXSCylinders](http://mcxtrace.org/download/components/3.1/contrib/SAXSCylinders.html), [SAXSPDB](http://mcxtrace.org/download/components/3.1/contrib/SAXSPDB.html) and many more for SAXS.
- [Powder](https://www.mcxtrace.org/download/components/3.1/samples/PowderN.html): diffraction.
- [Polycrystal](https://www.mcxtrace.org/download/components/3.1/samples/Polycrystal.html): diffraction.
- [Single crystal](https://www.mcxtrace.org/download/components/3.1/samples/Single_crystal.html): diffraction, also for MX.
- [Pump-probe](https://www.mcxtrace.org/download/components/3.1/samples/Molecule_2state.html) (2 states) molecule to simulate a laser-probe decay (time resolved).
- [Absorption_sample](https://www.mcxtrace.org/download/components/3.1/samples/Absorption_sample.html) a 1 or 2 absorbing materials as a box or cylinder; [Filter](https://www.mcxtrace.org/download/components/3.1/optics/Filter.html) which can handle absorption and refraction, as a block or any geometry; [Abs_objects](https://www.mcxtrace.org/download/components/3.1/samples/Abs_objects.html) a set of absorbing objects which geometry is set from OFF/PLY files.
- [Isotropic_Sqw](https://mcxtrace.org/download/components/3.1/samples/Isotropic_Sqw.html) an coherent inelastic scattering sample for isotropic density materials (e.g. amorphous, liquids, powders) that models Thompson scattering from a dynamic structure factor S(q,w) (IXS). Handles both elastic and inelastic contributions. Handles OFF/PLY anyshape geometry. This is an *experimental* sample.
- [Fluorescence](https://github.com/McStasMcXtrace/McCode/blob/mccode-3/mcxtrace-comps/samples/Fluorescence.comp) handles absorption, fluorescence, Compton and Rayleigh scattering, for any chemical formulae. Handles OFF/PLY anyshape geometry. This is an *experimental* sample (using [XRayLib](https://github.com/tschoonj/xraylib/wiki)).

All samples can have simple geometric shapes (some incl. hollow shapes).
Powder, SX, Filter, Abs_objects, Sqw and Fuorescence can have any shape (PLY/OFF).
Powder sample supports multiple concentric geometries (e.g. for cryostat, containers, ...).
McXtrace comes with a [material data base](https://www.mcxtrace.org/download/components/3.1/data), and can use e.g. NIST files.

## Monitors

There are plenty of monitors. Some are specific (e.g. a [PSD](https://www.mcxtrace.org/download/components/3.1/monitors/PSD_monitor.html), others are versatile (e.g. [Monitor_nD](https://www.mcxtrace.org/download/components/3.1/monitors/Monitor_nD.html)).

## Exercise: look at a component and modify it

In this hands-on, we shall copy an existing component, and change its behaviour. We start from the Session 1 *[First beam-line](../2_1st_Beamline)*. The monitor used is *[PSD_monitor](https://www.mcxtrace.org/download/components/3.1/monitors/PSD_monitor.html)*. 

This is a very simple example. When you get more experienced with McXtrace, we encourage you to contribute by sending to us (via [GitHub](https://github.com/McStasMcXtrace/McCode) or email at mcxtrace-users@mcxtrace.org) your own components, so that they become available to others.

Start by copying the initial component to your working directory, and edit it:
```bash
cp /usr/share/mcxtrace/3.1/monitors/PSD_monitor.comp myPSD_monitor.comp
gedit myPSD_monitor.comp &
```

### Step 1: A component dissection

We find the same structure as a beam-line description, e.g. `DEFINE`, `INITIALIZE`, `TRACE`.

The `DEFINE` section of the component contains its visible part, that is its name, and the parameters that control its behaviour. If this case, we mostly find dimensioning and binning parameters, with their default values.

The `INITIALIZE` section is executed only once at the start of the simulation. In this rather complex example, we check the geometry (rectangle, circle, etc...), define internal variables and allocate arrays.

The `TRACE` section is where we actually *do* something with the rays. An initial ray is described with variables:
- `x,y,z, kx,ky,kz, phi, t, Ex,Ey,Ez, p`

In the `PSD_monitor` component, there is a test for intersection of the detector surface (remember axis *X* is e.g. left, *Y* is vertical, *Z* is forward). When intersection is found, the statistical weight *p* is added to the monitor intensity `PSD_p`. The `PSD_N` holds the number of recorded events (per bin), and `PSD_p2` holds the square of the statistical weight, used to compute the error bars on intensity.

The `SAVE` section is where the output file is generated. In most cases, the `DETECTOR_OUT_2D` function is called with name, labels, bounds, arrays, and file name. This optional section is mostly used in monitor components.

The `FINALLY` section is used here to free the allocated memory on exit, and the `MCDISPLAY` indicates how to draw the component with drawing primitives (line, multiline, rectangle, box, circle, sphere, off_display).

### Step 2: Component surgery

Let's now come back to the `TRACE` section. In this exercise, we wish to make the monitor energy sensitive. A paper from James Holton (https://bl831.als.lbl.gov/~jamesh/mlfsom/Holton_R-factor_gap_2014_supp.pdf)  shows that in first approximation the number of visible photons generated on a Tb:Gd<sub>2</sub>O<sub>2</sub>S phosphor layer is about E<sub>ph</sub>/25 where E is the photon energy in eV. This means that statistically, an incident X-ray photon decays into 25 eV visible photons.

Identify the variable `e` and use it to compute the number of visible photons using the detected photon energy. *Tip*: E<sub>ph</sub> = K2E |k|, where K2E is a defined conversion factor within McXtrace. The resulting energy is in keV.

Modify the statements (for rectangle and circle cases) where the weight `p` is accumulated, and use `e` instead.

### Step 3: Simulating a more realistic detector response

Make a copy of the *First beam-line* description that was assembled in the previous hands-on. Change the beam-line description to use the newly created component, and run it.

Does it change the overall response ? The effect can mostly be seen when changing the photon energy. To achieve that, the source enegry must be made tunable. 

Change the `DEFINE ... (Par1=1 )` line into DEFINE ... `(E0=15)`. We turn the `Par1` input parameter into something meaningful. Change the Source component in the TRACE section to use `E0` as central photon energy (in keV). Now run the simulation again, and request a scan of the `E0` parameter. For this you need to specify in the Run dialogue a scanning range as `E0=1,50` and a number of steps (_Sweeps_), say 25.

Plot the results and comment.

:bulb: You may as well emit a white/rose spectrum at the source (with a large `dE`) and use a similarly modified `E_monitor` to acquire in a single computation the whole energy range and detector response. Using an unmodified `E_monitor` can be used to compare the ideal detector output with your version.

---
*McXtrace training - 2023*
