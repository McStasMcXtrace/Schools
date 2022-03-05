# McXtrace training: samples and virtual experiment: tomography

In this practical session we shall simulate a very simple model of an absorption spcetroscopy beam-line. We shall use a sample with an "any shape" complex volume, which we can rotate to simulate a sinogram. We shall also simulate a (very small) hyper-spectral data set.

## Absorption data files
There is a dedicated documentation and tool to get absorption data files. 
- [HOWTO: McXtrace absorption files (materials)](https://github.com/McStasMcXtrace/McCode/wiki/HOWTO%3A-McXtrace-absorption-files-%28materials%29)

Usual materials are already available in the [data](http://mcxtrace.org/download/components/3.0/data/) directory.
:warning: However, we currently only treat monoatomic elements, without handling the structure (i.e. no EXAFS yet).

## A simple absorption/tomography station

There are a few absorption sample components that can be used:
- [Absorption_sample](http://www.mcxtrace.org/download/components/3.0/samples/Absorption_sample.html) a 1 or 2 absorbing materials as a box or cylinder.
- [Abs_objects](http://www.mcxtrace.org/download/components/3.0/samples/Abs_objects.html) a set of absorbing objects which geometry is set from OFF/PLY files.
- [Filter](http://www.mcxtrace.org/download/components/3.0/optics/Filter.html) which can handle absorption and refraction, as a block or any OFF/PLY geometry.

A typical beam-line should look like:
- a photon source
- some optics/slits to shape the beam
- a rotating stage carrying a sample
- a detector

A good start is to search for examples that already use an absorption sample. We find the [`Airport_scannerII.instr`](http://mcxtrace.org/download/components/3.0/examples/Airport_scannerII.html) in group DTU, and the [`SOLEIL_ROCK`](http://mcxtrace.org/download/components/3.0/examples/SOLEIL_ROCK.html) beam-line. 

Lets start with the [`Airport_scannerII.instr`](http://mcxtrace.org/download/components/3.0/examples/Airport_scannerII.html) in which we substitute the `sample_scan` component with a `Filter` one (the `Abs_objects` seems broken). Let's use a geometry file `wire.ply` made of Mn.

:runner: Start a computation with 1e6 photon events, possibly with MPI. Plot it.

<img src="images/Airport_scanner.png">

# tomogram
 
:runner: Now, do a rotation of the sample around the vertical axis with `ANGLE=0,180` in 10 steps. Use 1e6 photon events, and MPI.

# absorption edge

To make sure absorption is properly taken into account, we shall now perform an energy scan. Looking at the current model, the energy spread is white. Lets restrict it to `E0` (default 6 keV), and `dE=0.1`. For this, a input parameters `E0=6` and `dE=0.1` should be set in the DEFINE INSTRUMENT line. We shall scan `E0`. This way we can change the white beam into a monochromatic distribution easily.

Also, the beam focus should be transformed into a small area. **TO BE DONE...**

:runner: Perform an energy scan between 5.8 and 6.8 keV. Plot the results.

# hyper-spectral imaging

With a narrow beam, we shall record the XANES as a function of the position of the beam hitting the sample. For this, we shall scan the sample position, and measure the outgoing spectrum for a flat incoming spectrum. To get a faster simulation, we use a white beam around 6 keV.

:runner: Add an energy monitor right after `psd2`.

Happy simulating!

---



---
*McXtrace training - 2022*
