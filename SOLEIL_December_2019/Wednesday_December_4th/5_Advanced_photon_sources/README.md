# Source modelling.
In this pratical exercise we will try out a few ways of modelling an X-ray source in McXtrace and some possibilitiers for coupling McXtrace to other packages.

## Exercise: Using the native McXtrace Undulator model.
N.b. This requires the Gnu Scientific Library to be installed. This should have been automatically installed with McXtrace install procedure. Should for some reason this not be the case you may go to (https://github.com/McStasMcXtrace/Schools/wiki/GSL-Installation) to find installation instructions.
1. Start a new simulation and insert an Undulator source in it. The Undulator component has many possible parameters. To be able to compare with later results we'll use what corresponds to the default setting of the s.k. "standard undulator" at SPring-8 in Japan: <code>
    E0=13, dE=1, Ee=8, dEe=0.001, Ie=0.1, 
    K=1.03118, Nper=140, lu=3.2e-2, 
    sigey=6.17e-6, sigex=0.29979e-3, sigepx=0.01226e-3, sigepy=1.1e-6, 
    focus_xw=1e-4, focus_yh=1e-4, dist=20, 
    E1st=12.400
</code>
Just as you have seen  before, the focus_xw, focus_yh, and dist parameteres simply indicate at sampling window downstream of the undulator. In McXtrace, the (0,0,0)-point is taken to be the exit plane of the undulator.
2. Insert two monitors **20 m** downstream: one PSD and one energy-resolved monitor. Make sure that the monitors are big enough to catch all the radiation you expect, including the energy range. 

**Unfortunately there's a 1-off bug in Source_spectra.com** Please use this ![.comp file](data/Source_spectra.comp?raw=true)

## Exercise: Connect with SPECTRA
If you do not have it already you may download SPECTRA freely from the riken website: (http://spectrax.org/spectra/), but for the purpose of this exercise we have pre-generated a set of datafiles, that you may use: [1st harmonic](data/sp8sU_h1.zip?raw=true ""), and [3rd harmonic](data/sp8sU_h3.zip?raw=true "").

### Details for the actual data files:
1. Generated using SPECTRA 10.0, using the Standard Spring-8 linear Undulator model, and standard beamline settings, I.e. as generic as possible. The opening screen look like this:
![spectra mani screen](images/spectra_main.png?raw=true "")

2. The type of calculation that McXtrace expects is:`Photon distribution at source point -> Wigner function -> Phase-space profile -> x-x' Plane (Projected)`, and correspondingly for the `x-y'` plane. These two filesets may then be used to drive a source in McXtrace.
3. Energy ranges may be generated automatically using the detuning parameter. 

There are some limitations in this incarnation:
1. The source in McXtrace is modelled as a point source.
2. Only a single harmonic is included. Stitching could solve this.
3. The projected planes may not be representative of the off-peak field. 

### Set up a simulation of this.
Download and unzip the files. In both cases the filename encodes the relevant energy interval:
1st harmonic: Emin=7.44 keV Emax=15.8 keV
3rd harmonic: Emin=37 Emax=58 kev TODO check these values

* Start a new simulation .instr file and insert a Source_spectra. Set parameter values for 
**spectra_stem_x=""sp8sU_h1_e7p44_18p6_x/sp8sU_h1_e7p44_18p6_x"
 spectra_stem_y=""sp8sU_h1_e7p44_18p6_y/sp8sU_h1_e7p44_18p6_y"  
nx=21 npx=21 ny=21 npy=21, nE=11, Emin=7.44, Emax=18.6**. Also set **E0** and **dE** to something that fits within **Emin** and **Emax**
* Insert monitors downstream of the source to monitor the source radiation. Catch the radiation on an energy resolved monitor see what the peak looks like.
* Try to move a point-like energy resolved detector around the radiation field. Does it beahve as you thing.
* Insert a slit 5 m downstream of the source on the optical axis. This mimicks the 1st order behaviour of front-end apertures. Scan the pinhole size to investigate the energy spectrum as a function of slit opening.  
* Focus the beam again using one of the techniques from yesterdays exercise. 

## Exercise: Use SRW-generated output to simulate a beamline.

**N.B. The gremlins that be have introduced a bug in MCPL_input.comp - there is an update version in the code-tree and ![here](data/MCPL_input.comp)**

We will now use a different utiltiy to drive a McXtrace-simulation: MCPL. In this case the MCPL-file is actually generated using [SRW](https://www.github.com/ochubar/SRW). This may be done using a c++-program `srw2mcl` that makes repeated calls to SRW and generates rays from that. This procedure is rather slow, so for this tutorial (to save time) we provide a pre-generated mcpl-file that you may use. The program itself may be found in the McXtrace repository and on the  McXtrace website, and you need a working installation of SRW to compile it.

The file you need is called [sp8stdU.mcpl.gz](data/sp8stdU.mcpl.gz?raw=true ""). There is also a bigger version of this same file for better sampling [sp8stdUl.mcpl.gz](data/sp8stdUl.mcpl.gz?raw=true ""). But this obviously takes longer to download.

* Use the McXtrace component **MCPL_input** to read rays from it and start them in a McXtrace simulation. 
* Insert a downstream `PSD_monitor` and an `E_monitor` to catch the generated radiation. Leave some room (2 m or so) between the mcpl-file and the monitor. SRW considers the undulator centre its reference point and so rays may actually originate there. 
* What was the fundamental energy of the 1st harmonic?
* This procedure relies on the undulator spectrum being sufficiently sampled by the srw2mcpl-program. Determine the sampling limits of the file using your monitors. 
* Open your "old" instrument from before, and replace the source with the MCPL_input solution.
