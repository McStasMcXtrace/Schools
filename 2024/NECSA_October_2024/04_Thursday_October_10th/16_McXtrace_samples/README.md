# McXtrace Samples

![McXtrace](../../pics/mcxtrace-logo.png  "McXtrace")

Samples are essential to build so-called virtual beam-lines, and reproduce data that look like real experiments. There are dedicated sessions on this topic.

To have an overview of existing sample, have a lokk at the [official samples](https://github.com/McStasMcXtrace/McCode/tree/main/mcxtrace-comps/samples) and [contributed samples](https://github.com/McStasMcXtrace/McCode/tree/main/mcxtrace-comps/contrib).

## Overview

Here is a list of the main components you may use:

Small angle scattering (large molecules/structures, polymers, colloids, ...)

- [SaSView models](https://github.com/McStasMcXtrace/McCode/tree/main/mcxtrace-comps/sasmodels): 96 models from SasView, PDB, Nanodiscs, Liposomes, I(q), ... [Saxs_sphere](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Saxs_spheres.comp), [SAXSCylinders](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/contrib/SAXSCylinders.comp), [SAXSPDB](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/contrib/SAXSPDB.comp) and many more for SAXS.

Diffraction (ordered crystals, including proteins)

- [Powder](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/PowderN.comp): diffraction. Can read CIF files via `cif2hkl`. One of the most efficient and versatile sample.
- [Polycrystal](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Polycrystal.comp): diffraction (quite complex to use).
- [Single crystal](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Single_crystal.comp): diffraction, also for MX. Can read CIF files via `cif2hkl`. Very efficient, can also model powders and some textures.

Spectroscopy (scattering depends on energy)

- [Fluorescence](https://github.com/McStasMcXtrace/McCode/blob/mccode-3/mcxtrace-comps/samples/Fluorescence.comp) handles absorption, fluorescence, Compton and Rayleigh scattering, for any chemical formulae (incl CIF files). Handles OFF/PLY anyshape geometry. This component is based on [XRayLib](https://github.com/tschoonj/xraylib/wiki).
- [Pump-probe](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Molecule_2state.comp) (2 states) molecule to simulate a laser-probe decay (time resolved).
- [Absorption_sample](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Absorption_sample.comp) a 1 or 2 absorbing materials as a box or cylinder; [Filter](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/optics/Filter.comp) which can handle absorption and refraction, as a block or any geometry; [Abs_objects](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Abs_objects.comp) a set of absorbing objects which geometry is set from OFF/PLY files. The `Fluorescence` sample has the same capabilities for any material, and may be a very handy alternative.
- [Isotropic_Sqw](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Isotropic_Sqw.comp) an coherent inelastic scattering sample for isotropic density materials (e.g. amorphous, liquids, powders) that models Thompson scattering from a dynamic structure factor S(q,w) (IXS). Handles both elastic and inelastic contributions. Handles OFF/PLY any shape geometry. This is an *experimental* sample.

All samples can have simple geometric shapes (some incl. hollow shapes).
`PowderN`, `Single_crystal`, `Filter`, `Abs_objects`, `Isotropic_Sqw` and `Fluorescence` can have any shape (PLY/OFF).
Powder sample supports multiple concentric geometries (e.g. for cryostat, containers, ...).
McXtrace comes with a [material data base](https://github.com/McStasMcXtrace/McCode/tree/main/mcxtrace-comps/data), and can use e.g. NIST files.

It is best to first search for existing beam-line models that make use of these samples, to learn how to configure and insert the sample in a model.

## Absorption (spectroscopy and tomography)

The absorption spectroscopy is a very simple measurement technique. The idea is to send an X-ray beam (white, pink or monochromatic beam), and illuminate a sample. The incident X-ray photons then traverse the sample volume. The absorption fraction depends on the incident energy. Indeed, above a given energy for each atom (the threshold), the X-rays eject inner electrons (e.g. from the K-edge, photo-emission). The energy levels are perfectly tabulated and specific to each atom. The X-rays are then 'absorbed' which means that the transmitted beam is decreased.
The transmission follows the classical Beer's law:

$I/I_0 = exp^{-d \mu(E)}$

where $I_0$ is the incoming intensity, $I$ is the transmitted intensity, $d$ is the propagation distance into the material, and $\mu(E)$ is the absorption coefficient. In practice, the data is normalised, and we rather show 1-T.

So, by just changing the incident energy, across absorption edges, it is possible to identify the material composition (XAS), as well as its oxidation state (XANES) and even local neighbours (EXAFS). This is the absorption spectroscopy.

Illuminating a volumetric sample, and placing an image detector after the sample, a transmitted projection is obtained. The images are 'semi-transparent' as a function of the X-ray energy and material. By rotating the sample, and taking many images, it is possible to reconstruct the 3D volume (with the object internals) from the projections. At Synchrotron SOLEIL, we use codes such as PyHST2, Nabu, TomoPy, Astra, and UFO. This is the tomography.

The following image has been obtained with the `Test_Absorption` model, which goes through a block of manganese Mn. A polychromatic beam goes through the sample. The energy-sensitive detector shows the absorption spectra without the need to scan (one of the many advantages of simulations).

![Absorption Mn](pics/Absorption.png)

The top curve shows intensity as a function of the energy. There is a drop after the Mn K-adge $E_K=6539$ keV. The image bellow shows the shadow of the block. The absorbed X-rays are converted into e.g. fluorescence and Auger electrons (not modelled here).

## Fluorescence

![Fluo](pics/Fluo.png)

## Powder diffraction

![PowderN](pics/PowderN.png)

## Single crystal diffraction

![SX](pics/SX.png)

## Small angle-scattering (diffraction)

![SAXS](pics/SAXS.png)


