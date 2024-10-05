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
- [Absorption_sample](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Absorption_sample.comp) a 1 or 2 absorbing materials as a box or cylinder; [Filter](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/optics/Filter.comp) which can handle absorption and refraction, as a block or any geometry; [Abs_objects](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Abs_objects.comp) a set of absorbing objects which geometry is set from OFF/PLY files. The `Fluorescence` sample has the same capabilities for any material, and should be prefered still.
- [Isotropic_Sqw](https://github.com/McStasMcXtrace/McCode/blob/main/mcxtrace-comps/samples/Isotropic_Sqw.comp) an coherent inelastic scattering sample for isotropic density materials (e.g. amorphous, liquids, powders) that models Thompson scattering from a dynamic structure factor S(q,w) (IXS). Handles both elastic and inelastic contributions. Handles OFF/PLY anyshape geometry. This is an *experimental* sample.

All samples can have simple geometric shapes (some incl. hollow shapes).
`PowderN`, `Single_crystal`, `Filter`, `Abs_objects`, `Isotropic_Sqw` and `Fluorescence` can have any shape (PLY/OFF).
Powder sample supports multiple concentric geometries (e.g. for cryostat, containers, ...).
McXtrace comes with a [material data base](https://github.com/McStasMcXtrace/McCode/tree/main/mcxtrace-comps/data), and can use e.g. NIST files.

It is best to first search for existing beam-line models that make use of these samples, to learn how to configure and insert the sample in a model.
