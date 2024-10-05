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