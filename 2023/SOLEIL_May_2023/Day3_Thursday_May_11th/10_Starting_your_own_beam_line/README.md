# Starting your own beam-line model: we help you !

In this 'anarchy' session, we shall do our best to get you started with your own beam-line model.

## Global methodology

1. Assemble a succinct list of elements that compose your beam-line. As a start, focus on the main components.
2. Identify which [McXtrace components](http://mcxtrace.org/download/components/) best match your beam-line elements. Prefer simple components when there are many possibilities.
3. Search for existing models that use these components. They are often in the 'Tests' category, but you may as well find a model which is pretty close in spirit to your beam-line. These examples will demonstrate how to actually use the components (syntax, parameters, etc) you envisage. Get inspiration here.
4. For each identified [McXtrace component](http://mcxtrace.org/download/components/), look at their main properties (parameters, geometry) and search for these in real life. This means you will need to estimate plausible values, and potentially refer to some documentation. As a start, use approximate values to avoid loosing time in details.
5. Position these elements in space w.r.t. the other elements. :bulb: A special note on mirrors and monochromators: the positioning/orientation of these can be tricky. It is usually a good advice to search for examples that contain them, and copy/paste/adapt the corresponding code in your own model. You will save a lot of time scratching your head this way.
6. Once done, you may optionally identify which values/parameters could be tunable in the beam-line model in order for instance to further scan or optimize their values. This can be an initial energy, a slit aperture size, a sample tilt, etc. These can then be moved to the DEFINE INSTRUMENT line as a variable, and transferred where appropriate in the INITIALIZE and TRACE sections. Set default values for these parameters.
7. Complete the header at the beginning of the `.instr` file to reflect your institution (SOLEIL), the beam-line, yourself as author, a description of the beam-line (copy-paste from [SOLEIL website](https://www.synchrotron-soleil.fr/en/beamlines)), a description of the model parameters (those from the step above).
7. Send us your beam-line, even basic/incomplete, so that we can further help you and enrich the gallery of SOLEIL beam-lines given as examples. By gathering beam-line models and distributing them with McXtrace releases, we are all happy to find examples to start with. 

## Beam-line model structure

The structure of a beam line should have:

1. A photon source
2. Some optics (e.g. slits, monochromators, mirrors, lenses, ...)
3. A sample model as it provides a result which goes beyond a simple photon beam. Using even a simple sample component makes the model much more convincing scientifically, and it also allows to estimate the photon counts at the detector. The `PowderN` `Single-crystal`, `Saxs_spheres.comp` or `Absorption_sample` are good bets.
4. Monitors (you may position as many as required, not only at the end). The `Monitor_nD` is extremely versatile, but often you will start with e.g. `E_monitor` or `PSD_monitor`.

## Compilation errors

Of course, your model will hardly compile initially. Read the compilation output to decode the cryptic messages and identify where errors could be. An other strategy is to comment (use `/* ... */`) most of the components and leave only the source. Then compile and iteratively un-comment components one-by-one with intermediate compilation checks.

---

## SOLEIL beam-lines

For the SOLEIL photon sources, refer to the session 4 ["Sources"](https://github.com/McStasMcXtrace/Schools/tree/master/2023/SOLEIL_May_2023/Day2_Wednesday_May_10th/4_Sources) to get the electron beam parameters.

### SOLEIL LUCIA Fluorescence

A simplified LUCIA model would include the following elements:

![LUCIA layout](https://www.synchrotron-soleil.fr/sites/default/files/resize/lignes_de_lumieres/lucia/lucia_short-800x373.jpg)

Position | Element
---------|----------
0        | the HU52 undulator (see practical 4 "Sources") 2.75 GeV, 32x52.4 mm, e-beam 218.2x8 µm, 30x3.7 µrad
13.4     | 0.5 x 0.5 mm slit
31.5     | a DCM (Si111, theta=5-75 deg)
39-39.3  | a KB mirror set
39.63    | a sample stage with a Fluorescence component
39.7     | a set of detectors (e.g. XRF)

Reference:

- D. Vantelon et al, J. Sync Rad 23 (2016) 635; DOI: 10.1107/S1600577516000746
- https://www.researchgate.net/publication/294278058_The_LUCIA_beamline_at_SOLEIL

### SOLEIL DIFFABS diffraction/absorption

The [DIFFABS](https://www.synchrotron-soleil.fr/fr/lignes-de-lumiere/diffabs) beam-line at SOLEIL is illuminated with a Bender, such as the one discussed in session 4 "Sources".

![SOLEIL_DIFFABS](images/SOLEIL_DIFFABS.png)

In order to assemble a DIFFABS model, we may set-up the following elements:

Position | Element
---------|----------
0     | the BM D13-1 Bender B=1.71 T in range 3-23 keV, e-beam cross-section 55.1x20.6 µm2 123.7x1.6 µrad. Other specs give 60.1x24.9 µm and 134.8x2.1 µrad. Critical energy 8.6 keV.
11.85 / 12.67 | primary slit
14.92 / 14.83 / 14.922| M1 bent mirror with Rh/Si coating. Reduces the vertical divergence. length 1300 mm, width 100 mm. 2.4 mrad incidence above 19 keV, 3 mrad in 12.5-19 keV, and 5 mrad in 6.5-12.5 keV, and 6 mrad below.
17.462 | a Si(111) DCM. First crystal is 200x100 mm flat at 2.54 m from M1; 2nd is bent 70x100 mm to reduce the horizontal divergence (perp. to beam). See the practical 5 "Optics". 
19.28 | M2 bent mirror with Rh coating, 1300x100 mm. Focuses horizontally the beam, and increases the divergence. At 1.817 from DCM.
?     | _a KB mirror set for micro focusing. Optional, see the practical 5 "Optics"._
?     | _a Fresnel zone plate. Optional, we ignore it here._
31.45 | sample stage with e.g. Fluorescence and/or PowderN components. At 13.99 from DCM.
32.09 | a set of detectors (transmission, XRD and XRF), radius 654 mm

Reference: 

- Gallard, Thèse, (2019) ["Etude in situ de la cristallisation et des contraintes dans des nanostructures de GeTe par diffraction du rayonnement X synchrotron"](https://www.theses.fr/2019AIXM0037.pdf)

### SOLEIL PSICHE beam-line tomography/diffraction

![PSICHE](https://github.com/McStasMcXtrace/Schools/raw/master/2023/SOLEIL_May_2023/Day3_Thursday_May_11th/9_Practical_Virtual_Exp_tomo/images/SOLEIL_PSICHE.png)

Position  | Element
----------|----------
0         | 15-100 keV Wiggler 2.1 T lu=50 mm 41 periods (K=10) e-beam 333.3x5.9 µm, 17x5.3 µrad
?         | primary mirror (c1/c2), ignored here
17.5 m    | removable DCM Si(111) for diffraction contrast tomography
?         | KB mirrors for a focused beam 100x100 µm (ignored here as we work in white beam)
21 m      | sample area, beam size 16.8x5.9 mm2
21.5      | a set of detectors (transmission, diffraction, fluo)

References:

- https://www.synchrotron-soleil.fr/en/beamlines/psiche
- A. King et al, Rev Sci Instrum 87, 093704 (2016), DOI: 10.1063/1.4961365
- E. Boulard et al, J Sync rad 2018 25, 818, DOI: 10.1107/S1600577518004861
- Guignot et al 2013, https://ui.adsabs.harvard.edu/abs/2013AGUFMMR31A2287G/abstract

### SOLEIL SAMBA absorption

The photon source is the bending magnet D09-1 used in range 4-40 keV. The beam shape is defined by 2 focusing mirrors and a DCM.

Position  | Element
----------|----------
0         | bending magnet D09-1 e-beam 53.2x23.6 µm 119.4x1.5 µrad
12.8      | primary slits
14.1      | focusing cylindrically bent mirror aceptance=6.2 mrad 1200x88 mm, Pd coating, theta=1-10 mrad.
16.1      | DCM Si(220)
18.1      | collimating cylindrically bent mirror
30        | sample stage, beam 2x1 mm -> 300x200 um when focused with mirrors
?         | fluorescence detector at 90 deg, transmission in beam

Reference:

- V. Briois et al, UVX 2010 (2011) 41–47; DOI: 10.1051/uvx/2011006
- S Belin et al, Phys. Scr. 2005 980; DOI: 10.1238/Physica.Topical.115a00980


### SOLEIL SWING SAXS

The [SWING](https://www.synchrotron-soleil.fr/fr/lignes-de-lumiere/swing) beam-line is using a U20 undulator between. 

![SWING layout](https://github.com/McStasMcXtrace/Schools/blob/master/2023/SOLEIL_May_2023/Day3_Thursday_May_11th/9_Practical_Virtual_Exp_tomo/images/SWING.png).

The optics are basically:

Position | Element
---------|----------
0     | a U20 undulator (used between 5 and 16 keV) with sigma=321 (H) x 9.4 (V) µm and divergence 17 (H) x 3.5 (V) µrad
11.7  | a diaphgram (slit) 1x0.5mm2
20    | a Si (111) double monochromator
22.5  | a KB mirror set
31    | a CRL (f=81 cm)
32    | the sample position
+0.5-6.5 | a 162x155 mm EigerX4M detector, from sample

Then make use of the KB example from e.g. the [test_KB](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/Test_KB.instr). 
Use theta = 3 mrad and compute the KB curvatures to match the distance to the sample as f = R sin(theta/2). 

Then introduce the DCM from the [Template_DCM](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/Template_DCM.instr) (satisfying `TTH=RAD2DEG*asin(12398.42/(2*DM*E0))` with DM=5.4309 for Si111)

We suggest that, to start with, you skip the CRL, but if you so wish, you may extract the relevant CRL description from [Test_CRL](https://raw.githubusercontent.com/McStasMcXtrace/McCode/master/mcxtrace-comps/examples/Test_CRL_Be.instr). Adapt the number of lenses to get a focusing at about 1 m.

Add a sample as above, and a detector (such as `Monitor_nD(options="x y", bins=2000)`).

You should then get a simplified SWING model.


