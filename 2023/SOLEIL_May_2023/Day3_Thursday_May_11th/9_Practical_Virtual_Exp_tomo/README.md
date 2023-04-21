# McXtrace training: samples and virtual experiments: SAXS and tomography

In this session we shall simulate the output of simple models for:
- Small Angle X-ray Scattering (SAXS)
- Tomography (which is based on absorption), and diffraction

## Table of Contents
1. [Exercise A: Small angle scattering (SAXS)](#exercise-a-small-angle-scattering-saxs)
1. [Exercise B: Tomography](#exercise-b-tomograohy)

---

## Sample geometry (reminder from Samples/diffraction)

Sample components should be given a geometrical shape. The sample coordinate frame is usually (when not rotated) *X* on the left, *Y* vertical, and *Z* is 'forward'.

The geometry can be specified as:
- a sphere `radius=<value>`
- a cylinder `radius=<value>, yheight=<value>`
- box `xwidth=<value>, yheight=<value>, zdepth=<value>`
- any shape defined with a `geometry=<file>` with a [PLY](http://en.wikipedia.org/wiki/PLY_%28file_format%29)/[OFF](http://www.geomview.org/docs/html/OFF.html) file (vertices and polygons similar to STL). We provide example geometry files in the [data](http://mcxtrace.org/download/components/3.1/data/) directory (e.g. locally at `/usr/share/mcxtrace/x.y/data`). You may also use e.g. [Meshlab](https://www.meshlab.net/) or other geometry editors/modellers to create such files (rather simple text format). Not all samples support this geometry.

Some samples can be made hollow by specifying a `thickness` parameter. This is especially useful for containers (e.g. capillary) and sample environments.

Moreover, some samples support a `concentric` mode, which allows to insert a component within an other. We shall not consider this topic during this session.

Using the ROTATED keyword, you may orient this geometry in any direction.

---

## Exercise A: Small angle scattering (SAXS)

The small-angle X-ray scattering beam-lines measure very small beam deviations around the incident direction, following the Bragg-law _n_&lambda; = 2 _d_ sin(&theta;) where the incident wavelength &lambda; is fixed and we can see that small angles &theta; corresponds with large typical scattering unit sizes _d_.

<img src="https://www.researchgate.net/profile/Jingpeng-Li-2/publication/322112043/figure/fig1/AS:657582976954368@1533791407758/Schematic-diagram-of-the-incident-beam-of-SAXS.png">

There is large variety of SAXS sample models. Most of them correspond with isotropic scattering units.

The most complete one is using [SaSView models](https://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/models/index.html) from which about [60 have been ported](http://mcxtrace.org/download/components/3.1/samples/SasView_model.html) into McXtrace. These include isotropic and anisotropic models.

In the following, we shall start from the [TestSAXS](http://mcxtrace.org/download/components/3.1/examples/TestSAXS.html) example instrument (from the _Tests_). It models a toy SAXS beam-line with a set of possible sample models via the input parameter `SAMPLE` (some sample models have been inactivated or are buggy). It also has a PSD and a |Q| monitor (with radial integration).

To run this model, you will need a PDB files accessible at https://www.rcsb.org/structure/6LYZ. Connect to the site and retrieve the PDB file (top right -> Download File -> PDB).

#### Step A.1: simulate the scattering from a set of samples

Load the [TestSAXS](http://mcxtrace.org/download/components/3.1/examples/TestSAXS.html) beam-line model and open the 3D view (run in Trace mode).
Accumulate the photon rays (click on 'keep rays') and start to visualize the scattering pattern.

<img src="images/TestSAXS.png">

Now re-run in Simulation mode, with SAMPLE=0, 1, 4, and 11. Use MPI (recompile) with e.g. 4 cores and 1e6 rays. 
These correspond with:
- 0=SAXSSpheres
- 1=SAXSShells
- 4=SAXSLiposomes
- 11=SAXSPDBFast (can use PDB files to compute I(q))

Plot the results, and visualize the scattering curve of all samples.

<img src="images/SAMPLE_0.png" width="200" title="SAXSSpheres"> <img src="images/SAMPLE_1.png" width="200" title="SAXSShells">
<img src="images/SAMPLE_4.png" width="200" title="SAXSLiposomes"> <img src="images/SAMPLE_11.png" width="200" title="SAXSPDBFast">

:question: what can you say about the scattering units in the sample ? Do they compare/differ ?

#### Step A.2: simulate more complex samples

We now use the [TemplateSasView](http://www.mcxtrace.org/download/components/3.1/examples/templateSasView.html) in _Templates_.

As can be see, the default model index is number 10.

:question: 
- Identify which structure is being used by looking at the table [SasView_model](http://mcxtrace.org/download/components/3.1/samples/SasView_model.html). :warning: links are broken. You should follow the [SasView documentation](https://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/models/index.html). 
- Are we using the default SasView Cylinder model parameters ?

🏃 Run the simulation and plot the results.

A |q| detector would probably be a good idea. Add and instance of the `SAXSQMonitor` at 3.13 m away from the sample, with its `RadiusDetector=0.3`, the `DistanceFromSample=3.13`, `LambdaMin` and `Lambda0` set the nominal wavelength of the source, i.e. `lambda`.

🏃 Run the simulation again and plot the results, in Log-scale.

<img src="images/templateSasView.png" title="SasView_model cylinder">

Now change the structure to a bcc-paracrystal, using the default parameter values extracted from the [SasView documentation](https://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/models/index.html).

:runner: Run the simulation with the [bcc_paracrystal](https://www.sasview.org/docs/user/models/bcc_paracrystal.html) which is `SasView_model(index=4)`.

#### Step A.3: (optional) simplified SWING beam-line at SOLEIL

The [SWING](https://www.synchrotron-soleil.fr/fr/lignes-de-lumiere/swing) beam-line is using a U20 undulator between. 
The BL layout is available [here](images/SWING.pdf).

The optics are basically:

- a U20 undulator (used between 5 and 16 keV) with sigma=388 (H) x 8.1 (V) um and divergence 14.5 (H) x 4.6 (V) urad
- a diaphgram (slit) 1x0.5mm2 at 11.7 m from the source
- a Si (111) double monochromator at 20m from the source
- a KB mirror set at 22.5m from the source
- a CRL (f=81 cm) at 31m from the source
- the sample position at 32 m from the source
- a 162x155 mm EigerX4M detector at 0.5-6.5m from the sample

We suggest you first insert a U20 undulator, as described in the [SOLEIL_PX2a](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/SOLEIL_PX2a.instr) beam-line:
``` c
Undulator(
    E0=E0, 
    dE=1, 
    Ee=2.75, 
    dEe=0.001, 
    Ie=0.5, 
    K=1.788, 
    Nper=80, 
    lu=2.4e-2, 
    sigey=9.3e-6, 
    sigex=215.7e-6, 
    sigepx=29.3e-6, 
    sigepy=4.2e-6, 
    dist=29.5, 
    E1st=12.400)
```

Then make use of the KB example from e.g. the [test_KB](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/Test_KB.instr). 
Use theta = 3 mrad and compute the KB curvatures to match the distance to the sample as f = R sin(theta/2). 
``` c
COMPONENT mirror_curved = Mirror_curved(
    radius=R,
    length=1,
    width=1)
AT (0, 0, 31.5) RELATIVE PREVIOUS
ROTATED (0, RAD2DEG*theta, 0) RELATIVE PREVIOUS
EXTEND
%{ 
	if (!SCATTERED) ABSORB; 
%}

COMPONENT arm = Arm()
AT (0, 0, 0) RELATIVE PREVIOUS
ROTATED (0, RAD2DEG*theta, 0) RELATIVE PREVIOUS

COMPONENT arm_2 = Arm()
AT (0, 0, 1.5) RELATIVE PREVIOUS
ROTATED (0, 0, 90) RELATIVE PREVIOUS

COMPONENT mirror_2 = Mirror_curved(
    radius=R,   
    length=1,
    width=1)
AT (0, 0, 0) RELATIVE PREVIOUS
ROTATED (0, RAD2DEG*theta, 0) RELATIVE PREVIOUS
EXTEND
%{ 
	if (!SCATTERED) ABSORB; 
%}

COMPONENT arm_3 = Arm()
AT (0, 0, 0) RELATIVE PREVIOUS
ROTATED (0, RAD2DEG*theta, 0) RELATIVE PREVIOUS
```

Then introduce the DCM from the [Template_DCM](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/Template_DCM.instr) (satisfying `TTH=RAD2DEG*asin(12398.42/(2*DM*E0))` with DM=5.4309 for Si111)
``` c
COMPONENT dcm_xtal0 = Bragg_crystal(
    length=0.04, width=0.04, 
    alpha=alpha, h=1, k=1, l=1, material="Si.txt", V=160.1826)
AT(0,0,0.02) RELATIVE PREVIOUS
ROTATED (-TTH,0,0) RELATIVE PREVIOUS

COMPONENT dcm0 = Arm()
AT(0,0,0) RELATIVE dcm_xtal0
ROTATED (-TTH,0,0) RELATIVE PREVIOUS

COMPONENT dcm_xtal1 = COPY(dcm_xtal0)
AT(0,0,dcm_gap) RELATIVE dcm0
ROTATED (TTH,0,0) RELATIVE dcm0

COMPONENT dcm1 =Arm()
AT(0,0,0) RELATIVE dcm_xtal1
ROTATED (TTH,0,0) RELATIVE dcm_xtal1 
```

We suggest that, to start with, you skip the CRL, but if you so wish, you may extract the relevant CRL description from [Test_CRL](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/Test_CRL_Be.instr). Adapt the number of lenses to get a focusing at about 1 m.
``` c
COMPONENT lens_parab = Lens_parab(
    material_datafile = "Be.txt",
    r=200e-6, 
    r_ap=0.5e-3, 
    d=50e-6, 
    N=16)
AT (0, 0, 1) RELATIVE PREVIOUS
```

Add a sample as above, and a detector (such as `Monitor_nD(options="x y", bins=2000)`).

You should then get a simplified SWING model.


---


## Exercise B: Tomography

In this practical session we shall simulate a very simple model of an absorption spectroscopy beam-line. We shall use a sample with an "any shape" complex volume, which we can rotate to simulate a sinogram. We shall also simulate a (very small) hyper-spectral data set.

### Heating up: Absorption data files
There is a dedicated documentation and tool to get absorption data files. 
- [HOWTO: McXtrace absorption files (materials)](https://github.com/McStasMcXtrace/McCode/wiki/HOWTO%3A-McXtrace-absorption-files-%28materials%29)

Usual materials are already available in the [data](http://mcxtrace.org/download/components/3.1/data/) directory.
:warning: However, we currently only treat monoatomic elements, without handling the structure (i.e. no EXAFS).

### A simple absorption/tomography station

There are a few absorption sample components that can be used:
- [Absorption_sample](http://www.mcxtrace.org/download/components/3.1/samples/Absorption_sample.html) a 1 or 2 absorbing materials as a box or cylinder.
- [Abs_objects](http://www.mcxtrace.org/download/components/3.1/samples/Abs_objects.html) a set of absorbing objects which geometry is set from OFF/PLY files.
- [Filter](http://www.mcxtrace.org/download/components/3.1/optics/Filter.html) which can handle absorption and refraction, as a block or any OFF/PLY geometry.

A typical beam-line should look like:
- a photon source
- some optics/slits to shape the beam
- a rotating stage carrying a sample
- a detector

A good start is to search for examples that already use an absorption sample. We find the [`Airport_scannerII.instr`](http://mcxtrace.org/download/components/3.1/examples/Airport_scannerII.html) in group DTU, the [`NBI_Lab_TOMO`](http://mcxtrace.org/download/components/3.1/examples/NBI_Lab_TOMO.html) and the [`SOLEIL_ROCK`](http://mcxtrace.org/download/components/3.1/examples/SOLEIL_ROCK.html) beam-line. 

Lets start with the [`Airport_scannerII.instr`](http://mcxtrace.org/download/components/3.1/examples/Airport_scannerII.html) in which we substitute the `sample_scan` component with a `Filter` one (the `Abs_objects` seems broken). Let's use a geometry file `wire.ply` made of Mn. 


| :warning: If you use McXtrace 3.1, please copy the [3.x/Filter.comp](3.x/Filter.comp) fixed component |
|----|

:runner: Start a computation with 1e6 photon events, possibly with MPI. Plot it.

<img src="images/Airport_scanner.png">

### tomogram
 
:runner: Now, do a rotation of the sample around the vertical axis with `ANGLE=0,180` in 10 steps. Use 1e6 photon events, and MPI. Computation should last e.g. 1-3 minutes. Plot the results.

:runner: To visualize the individual images, use Ctrl-click on the `psd2_I` monitor. 

<img src="images/Airport_scanner-rotation.png">



### tomo/diffraction



----


