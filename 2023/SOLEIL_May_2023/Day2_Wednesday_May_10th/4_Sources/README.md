# Source modelling.
In this practical exercise we will try out a few ways of modelling an X-ray source in McXtrace and some possibilities for coupling McXtrace to other packages.

## Exercise: Using the native McXtrace Undulator model.

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/9f/Undulator.png/600px-Undulator.png">

1. Start a new simulation and insert an Undulator source in it. The Undulator component has many possible parameters. To be able to compare with later results we'll use what corresponds to the default setting of the s.k. "standard undulator" at SPring-8 in Japan: 
    <code>
    E0=13, dE=1, Ee=8, dEe=0.001, Ie=0.1, 
    K=1.03118, Nper=140, lu=3.2e-2, 
    sigey=6.17e-6, sigex=0.29979e-3, sigepx=0.01226e-3, sigepy=1.1e-6, 
    focus_xw=1e-4, focus_yh=1e-4, dist=20, 
    E1st=12.400</code>
    
Just as you have seen  before, the focus_xw, focus_yh, and dist parameters simply indicate at sampling window downstream of the undulator. In McXtrace, the (0,0,0)-point is taken to be the exit plane of the undulator.

2. Insert two monitors **20 m** downstream: one PSD and one energy-resolved monitor. Make sure that the monitors are big enough to catch all the radiation you expect, including the energy range.

There are also McXtrace components to model a `Wiggler` and a `Bending_magnet`. Last but not least, the [`Shadow_input`](http://mcxtrace.org/download/components/3.0/misc/Shadow_input.html) and [`Shadow_output`](http://mcxtrace.org/download/components/3.0/misc/Shadow_output.html) allow to read data events from [Shadow](https://github.com/oasys-kit/shadow3), but also to fully insert a McXtrace

:warning: N.B. This component requires the Gnu Scientific Library (GSL) to be installed. This should have been automatically installed with McXtrace install procedure. Should for some reason this not be the case you may go to (https://github.com/McStasMcXtrace/Schools/wiki/GSL-Installation) to find installation instructions. simulation as into a Shadow simulation sequence (read and write).

## Exercise: Connect with SPECTRA
If you do not have it already you may download [SPECTRA](http://spectrax.org/spectra/) freely from the Riken website, but for the purpose of this exercise we have pre-generated a set of datafiles, that you may use: [1st harmonic](data/sp8sU_h1.zip?raw=true ""), and [3rd harmonic](data/sp8sU_h3.zip?raw=true "").

### Details for the actual data files:
1. Generated using SPECTRA 11, using the Standard Spring-8 linear Undulator model, and standard beamline settings, i.e. as generic as possible. To get help about parameters, select the _Help_ menu, _Open Reference Manual_ item. The configuration file for this example is available as [spring8.json](data/spring8.json?raw=true ""). The opening screen look like this:
![spectra main screen](images/spectra_main.png?raw=true "")

2. The type of calculation that McXtrace expects is:`Photon distribution at source point -> Wigner function -> Phase-space profile -> x-x' Plane (Projected)`, and correspondingly for the `x-y'` plane. These two file-sets may then be used to drive a source in McXtrace. Choose ADCII format. Then select the _Run_ menu, _Start calculation_ item.
![spectra main screen](images/spectra_calculation_choice_x.png?raw=true "")

3. Energy ranges may be generated automatically using the de-tuning parameter (e.g. 0). Right-click on 'Detunning', select 'Scan this parameter' with e.g. 
- Initial	-0.4
- Final	0.5
- Number of points 11

Then select the _Run_ menu, _Start calculation_ item.

The data files are generated both in JSON, and ASCII. 

There are some limitations in this incarnation:

1. The source in McXtrace is modelled as a point source.
2. Only a single harmonic is included. Stitching could solve this.
3. The projected planes may not be representative of the off-peak field. 

### Set up a simulation of this.
Download and unzip the files. In both cases the filename encodes the relevant energy interval:

- 1st harmonic: Emin=7.44 keV Emax=15.8 keV
- 3rd harmonic: Emin=37 Emax=58 keV

Start a new simulation .instr file

* Insert a Source_spectra. Set parameter values for <code>
spectra_stem_x="sp8sU_h1_e7p44_18p6_x/sp8sU_h1_e7p44_18p6_x",
 spectra_stem_y="sp8sU_h1_e7p44_18p6_y/sp8sU_h1_e7p44_18p6_y",
 nE=11, Emin=7.44, Emax=18.6</code>  
 Also set **E0** and **dE** to something that fits within **Emin** and **Emax**
* Insert monitors downstream of the source to monitor the source radiation. Catch the radiation on an energy resolved monitor and see what the peak looks like.

The reason that you only see radiation in one quadrant is that the spectra datafiles (to save space) only contain data in this quadrant. By setting the parameters **symmetricx=1**, and **symmetricy=1**, the radiation field is mirrored in the ZX-, and ZY-planes respectively.

* Try to move a point-like energy resolved detector around the radiation field. Does it behave as you think it should.
* Insert a slit 5 m downstream of the source on the optical axis. This mimics the 1st order behaviour of front-end apertures. Scan the pinhole size to investigate the energy spectrum as a function of slit opening.  
* It is apparent that these supplied data-files are far too coarsely sampled. For a really useful simulation it is necessary to create bigger datafiles, for instance such as [1st harmonic long](data/sp8sU_h1_3.zip?raw=true "").  

Similarly, there is an interface with the [SIMPLEX](https://spectrax.org/simplex/index.html) and [GENESIS 1.3](http://genesis.web.psi.ch/index.html) codes for XFEL's.

## Exercise: Use SRW-generated output to simulate a beamline.

We will now use a different utiltiy to drive a McXtrace-simulation: MCPL. 

[MCPL](https://mctools.github.io/mcpl/) is an interchange file format to communicate with e.g. GEANT4, PHITS, MCNP, and SRW.

In this case the MCPL-file is actually generated using [SRW](https://www.github.com/ochubar/SRW). The expert tool for this purpose, a C++-program `srw2mcpl` is in development status and based on SRWlib and MCPL. The tool makes repeated calls to SRW and generates rays from that, and is not in its current form in an end-user state. 

The method of calling SRW pr. ray is rather slow, so for this tutorial (to save time) we provide a pre-generated MCPL-file that you may use. We are working on solutions to make the program available to McXtrace users in an easy-to use form. 

The file you need is called [sp8stdU.mcpl.gz](data/sp8stdU.mcpl.gz?raw=true ""). There is also a bigger version of this same file for better sampling [sp8stdUl.mcpl.gz](data/sp8stdUl.mcpl.gz?raw=true ""). But this obviously takes longer to download.

* Use the McXtrace component **MCPL_input** to read rays from it and start them in a McXtrace simulation. 
* Insert a downstream `PSD_monitor` and an `E_monitor` to catch the generated radiation. Leave some room (2 m or so) between the MCPL-file and the monitor. SRW considers the undulator centre its reference point and so rays may actually originate there. 
* What was the fundamental energy of the 1st harmonic?
* This procedure relies on the undulator spectrum being sufficiently sampled by the `srw2mcpl`-program. Determine the sampling limits of the file using your monitors. 
* Open your "old" instrument from before, and replace the source with the MCPL_input solution.

## Exercise: SOLEIL Undulator HU52 "Apple-II"

The LUCIA beam-line at Synchrotron SOLEIL (SD03C) is illuminated with an Undulator HU52 "Apple II" type (NdFeB magnets), 32 periods, gap 15 – 150mm, variable linear polarization, left and right circular polarizations, operating on harmonics 3 to 21. The energy range is 0.6-8 keV. 

The HU52 main parameters are:

HU52 parameter | symbol/unit | value
------------------------|-----------|-----
Period | (mm) | 52.4
Nb of periods  | lambda_U | 30 + 2
Gap | mm | 15.5-150
Field remanence | B(T) | 1.26
Magnetic field Z, max | Bz (T) | 1.974*exp(-3.1754*Gap/lambda_U)
Magnetic field X, max | Bx (T) | 1.901*exp(-4.3387*Gap/lambda_U)
Function beta horizontal |	beta_x (m)	|1.4
Function beta vertical|	beta_z ( m)	|1.4
emittance horizontale | ex (pm.rad)	|82.0
betatron	coupling | % | 30.0
emittance vertical |	ez (pm.rad)	|24.6
dimension RMS horizontal|	sigma_x (µm)	|10.7
dimension RMS vertical	|sigma_z (µm)	|5.9	
divergence RMS horizontal	|sigma_x' (µrad)	|7.7
divergence RMS vertical	|sigma_z' (µrad)|	4.2

The storage ring parameters for SOLEIL are available at https://www.synchrotron-soleil.fr/en/research/sources-and-accelerators/electron-beam-parameters/transverse-size-electron-beam-source
The bunch parameters are available at https://www.synchrotron-soleil.fr/en/research/sources-and-accelerators/electron-beam-parameters/bunch-length

We can then derive the Undulator component parameters for LUCIA on HU52:
```c
Undulator(
  E0 = 8,
  dE = 0.1,
  Ee = 2.75,
  dEe = 0.001,
  Ie = 0.5,
  tbunch = 47e-12,
  B = 0.42, // for a 15.5 mm gap
  Nper = 32,
  lu = 52.4e-3,
  sigey = 10.7e-6,
  sigex = 5.9e-6,
  sigepx = 7.7e-6,
  sigepy = 4.2e-6) 
```

For the DEIMOS beam-line (also using an HU52), we have the electron beam parameters:
- sizes of 188×8.2 μm RMS (H×V)
- divergences of 25.5×6.0  μrad RMS (H×V)
- we use a 30 mm gap corresponding to a magnetic field of 0.32 T

The detector is located at L=20.465 m from the photon source. The photon beam divergence at the detector should be around 1.2 μrad RMS.

References:
- F. Briquez et al., Proceedings of FEL08, Gyeongju, Korea, https://accelconf.web.cern.ch/fel2008/papers/tupph015.pdf
- Thierry Moreno et al., J Sync rad 19 (2012) 179, https://journals.iucr.org/s/issues/2012/02/00/kt5033/index.html
- https://accelconf.web.cern.ch/ipac2013/talks/mozb102_talk.pdf
- T. Moreno et al., https://www.researchgate.net/publication/258548494_Undulator_emission_analysis_Comparison_between_measurements_and_simulations


