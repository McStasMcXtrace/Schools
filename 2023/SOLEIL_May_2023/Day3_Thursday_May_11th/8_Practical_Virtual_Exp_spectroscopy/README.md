# McXtrace training: samples and virtual experiments: spectroscopy

In this session we propose to use some spectroscopy sample models:

- Absorption (using `Filter`, `Absorption_sample`, `Abs_objects`, and `Fluorescence` components)
- Fluorescence, togeteher with absorption, Compton and Rayleigh scattering
- Inelastic scattering (`Isotropic_Sqw`)

For this session you should start with:

- [Absorption, the ROCK beam-line](8b_Spectroscopy)
- [Fluorescence (experimental)](8c_Spectroscopy)
- [Inelastic scattering, a simple example](8a_Spectroscopy)

## Hyper-spectral imaging

We shall now use a complex geometry made with a material, and enclose it in a box with an other material. For this exercise, we shall use a GROUP arrangement. The 1st fluorescence sample will be set with a `geometry` parameter, while a 2nd Fluorescence component will be using a simple box with an other material. 

:question: can you identify why this is not a perfect solution ? What type of interaction process is missing ? Is there a way to be closer to the reality ?





With a narrow beam, we shall record the XAS signal as a function of the position of the beam hitting the sample. This corresponds with an hyper-spectral data set (X,Y,E).

In this exercise, we wish to demonstrate the concept by recording an image through the sample, as a function of the energy. In the current model, the incoming beam is wider than the sample material block, and the absorption contrast is not visible compared with the direct transmitted intensity. We thus choose to add a PSD monitor with dimensions matching the sample.

:runner: Add a 3x3 cm2 PSD monitor after the energy monitor. Perform an energy scan between 4.5 and 7.5 keV in 20 steps, with `dE=0.1` keV.

You should obtain a series of images, one per incoming energy. This is an hyper-spectral data set which can be merged into a 3D array with e.g. NumPy or Matlab. Of course, in this exercise, all images are similar as the sample is homogeneous. Putting back the `Filter` component with a complex geometry should allow to record a more realistic hyper-spectral data set.


---


