# McXtrace training: samples and virtual experiments: spectroscopy

In this session we propose to use some spectroscopy sample models:

- Absorption (using `Filter`, `Absorption_sample`, `Abs_objects`, and `Fluorescence` components)
- Fluorescence, togeteher with absorption, Compton and Rayleigh scattering
- Inelastic scattering (`Isotropic_Sqw`)

## Absorption edge

In the absorption process, a photon excites a core electron in an atom (e.g. K n=1 level), e.g. to a upper unoccupied level. When the photon energy is larger than the photo-ionization, the electron is ejected and it can be analyzed as in XPS, PEEM, ARPES etc. The emptied initial core level is then filled by upper level electron states (n=2,3...), and this transition generates a fluorescence line between levels. So, at the absorption edge, photons are absorbed, the transmitted intensity is decreased, and some lower desexcitation transitions emit lower energy X-ray fluorescence lines. 

![Absorption levels](https://upload.wikimedia.org/wikipedia/commons/b/bd/XASEdges.svg)

There are tables (see [RuppWeb](https://www.ruppweb.org/Xray/elements.html)) describing the absorption energies for each atom, and the corresponding fluorescence lines.

Only the `Fluorescence` component allows both any-spahe, and any compound material. The other absorption samples require to prepare absorption data files, and only work with pure mono-atomic materials.

## Assembling an absorption beam-line model

We here define a new beam-line built for absorption spectroscopy

1. Start a new beam-line, and set its input parameters as `E0`, `dE`, `theta` for the sample rotation angle, and the sample material and geometry (as string, not numbers). For instance `DEFINE INSTRUMENT SOLEIL_PSICHE(E0=6.5, dE=1, theta=0, string sample_material="MnFeCr", string sample_geometry="")`. We use a set of atoms with close absorption energies, as seen in the [edge energy tables](https://www.ruppweb.org/Xray/elements.html).

2. Insert a Wiggler such as the PSICHE@SOLEIL one as photon source:
``` c
Wiggler(E0 = E0, dE = dE, phase = 0, randomphase = 1, Ee = 2.4, Ie = 0.5, B = 2.1, Nper=41, sigey=9.3e-6, sigex=215.7e-6, length=38*50e-2)
```

3. Insert an energy and PSD monitor at 10 m from the source.

4. Insert a simple slit (2x2 mm) at e.g. 20 m away.

5. Insert an `Arm` component, as sample holder, at 0.5 m from the slit.

6. Add a `Fluorescence` sample on that `Arm`, and rotate it by `theta` along its vertical axis `Y` in order to be able to perform a tomography scan. Define its `material=sample_material` and `geometry=sample_geometry` using beam-line input argument (those defined in the DEFINE line above). Use a simple sample plate volume, e.g. 3x3x0.5 mm3.

7. Add a PSD detector at e.g. 10 cm after the sample, relative to the sample holder so that it does not also rotate with `theta`. Add as well an energy monitor (to catch the fluorescence), rotated by 45 deg wrt the incoming position, and an other one in transmission.

:runner: start a computation with a large energy range (e.g. `E0=11, dE=10`), and compare the incoming and transmitted flux. You should see the absorption edges for the given materials.

## Absorption spectroscopy

To make sure absorption is properly taken into account, we shall perform an energy scan as done on a real absorption spectroscopy BL. 

:runner: Let's restrict the energy range around `E0` with `dE=0.1`, and scan the energy between 4.5 and 8 in 20 steps. Plot the results. 

:question: Do you recover the transmission spectra ? Can you explain the differences ?

## Hyper-spectral imaging

We shall now use a complex geometry made with a material, and enclose it in a box with an other material. For this exercise, we shall use a GROUP arrangement. The 1st fluorescence sample will be set with a `geometry` parameter, while a 2nd Fluorescence component will be using a simple box with an other material. 

:question: can you identify why this is not a perfect solution ? What type of interaction process is missing ? Is there a way to be closer to the reality ?





With a narrow beam, we shall record the XAS signal as a function of the position of the beam hitting the sample. This corresponds with an hyper-spectral data set (X,Y,E).

In this exercise, we wish to demonstrate the concept by recording an image through the sample, as a function of the energy. In the current model, the incoming beam is wider than the sample material block, and the absorption contrast is not visible compared with the direct transmitted intensity. We thus choose to add a PSD monitor with dimensions matching the sample.

:runner: Add a 3x3 cm2 PSD monitor after the energy monitor. Perform an energy scan between 4.5 and 7.5 keV in 20 steps, with `dE=0.1` keV.

You should obtain a series of images, one per incoming energy. This is an hyper-spectral data set which can be merged into a 3D array with e.g. NumPy or Matlab. Of course, in this exercise, all images are similar as the sample is homogeneous. Putting back the `Filter` component with a complex geometry should allow to record a more realistic hyper-spectral data set. However :-1: it currently has a bug which prevents recovering the absorption spectrum. This will be fixed for the next release.


---


