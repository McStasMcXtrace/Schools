# Samples

McXtrace currently supports a set of samples:

- powder diffraction
- single crystal diffraction
- absorption spectroscopy (edge only, no pre-edge nor EXAFS)
- small angle scattering (with many models from SASview)
- fluorescence
- inelastic scattering (currently only isotropic density materials 2D IXS, no RIXS)

In the future we plan to add:

- photo-emission inelastic scattering (approx. 2D=XPS and 4D=ARPES)
- inelastic scattering in single crystals (4D IXS, no RIXS)

## Using "any-shape" geometry

As already discussed for the diffraction and tomography samples, some components support the "any shape" geometry. The component usually takes a `geometry=file` argument, for wich the `file` contains vertices and polygons stored into a [PLY](http://en.wikipedia.org/wiki/PLY_%28file_format%29)/[OFF](http://www.geomview.org/docs/html/OFF.html). These files can be generated from e.g MeshLab.

The following components currently support the OFF/PLY any-shape:

- `PowderN` powder diffraction
- `Single_crystal` single crystal diffraction
- `Filter` absorption
- `Abs_objects` absorption
- `Absorption_sample` absorption
- `Fluorescence` absorption, fluorescence, Compton, Rayleigh
- `Isotropic_Sqw` inelastic Thomson (IXS)
- `Shape` empty shape just for rendering in 3D
- `Monitor_Sqw` an inelastic scattering monitor
- `Monitor_nD` multi-purpose monitor

and we provide some example OFF/PLY in the 'data' directory from the "Docs".

## Combining samples

Some components allow to model hollow shapes:

- `PowderN` powder diffraction
- `Fluorescence` absorption, fluorescence, Compton, Rayleigh
- `Isotropic_Sqw` inelastic Thomson (IXS)

as well as multiple shapes/materials:

- `Abs_objects` absorption (a list of OFF/PLY shapes and materials)
- `Absorption_sample` absorption with an external and internal shape/material
- `Polycrystal` a map of crystal orientations

## Using the grammar for multiple samples/processes

#### WHEN 

One easy way to combine sample components is to use `WHEN`. The condition should still be set from e.g. a random variable. For instance, if you wish to combine a diffraction sample with a fluorescence/absorption one, you may use:

``` c
DECLARE %{
  double sample_choice;
%}

TRACE
...

COMPONENT sample_stage = Arm()
AT ...
EXTEND %{
  sample_choice = ceil(rand01()*2); /* rand01 selects a random number in 0-1 */
%}

COMPONENT sample1 = Fluorescence(...)
WHEN (sample_choice == 1)
AT (0,0,0) RELATIVE sample_stage

COMPONENT sample2 = PowderN(...)
WHEN (sample_choice == 2)
AT (0,0,0) RELATIVE sample_stage
```

This assumes the scattering probability for each sample is equal. You may scale the scattering with some EXTEND block after each sample:
``` c
EXTEND %{
  p = p/2; /* half scattering for 2 stacked scattering processes */
%}
```

#### GROUP

You may as well use the simpler `GROUP` keyword if:

- components do not "take" the whole photon events
- components are exclusive (one in a list)

Then you would typically write:
``` c
COMPONENT sample1 = Fluorescence(...)
GROUP samples_group
AT (0,0,0) RELATIVE sample_stage

COMPONENT sample2 = PowderN(...)
GROUP samples_group
AT (0,0,0) RELATIVE sample_stage
```

The choice is made when a component actually "SCATTER" (do something). All components are tested until one interacts, skipping the rest of the group.

#### Logic when interacting

Every component that "does something" sets the `SCATTERED` flag. 0 means nothing took place (no geometry intersection, no scattering).

So, a simple way to "remove" the non-interacting beam is to use:
``` c
EXTEND %{
  if (!SCATTERED) ABSORB;
%}
```
which is a perfect beam-stop for samples, only leaving the "scattered" part.

But you may perform more complex filtering based on the SCATTERED flag, or any other variable, best in conjunction with `WHEN`. You will find implementations in the `examples` (Docs).

