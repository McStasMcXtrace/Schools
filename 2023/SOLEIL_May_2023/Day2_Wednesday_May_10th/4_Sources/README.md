# Source modelling.

In this practical exercise we will try out a few ways of modelling an X-ray source in McXtrace and some possibilities for coupling McXtrace to other packages.

## Exercise: Using the native McXtrace Undulator model for SOLEIL photon sources

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/9f/Undulator.png/600px-Undulator.png">

The default Undulator model is well suited to model a photon source [Ref: K.J. Kim, AIP, conf. proc., 184, 1989. doi:10.1063/1.38046](https://pubs.aip.org/aip/acp/article/184/1/565/788822/Characteristics-of-synchrotron-radiation). 

Currently, there exists a set of undulator uses in the McXtrace examples.

Undulators | parameters
-----------|-----------
Test\_BM Undulator | `Undulator(Ee=1.5, K=1, E0=0.39, dE=0.2, Ie=0.4, B=0, gap=4.2, Nper=134, lu=3.65e-2, sigex=0.05367e-3, sigey=0.004e-3, focus_xw=10e-3,focus_yh=10e-3, dist=20)`
Test\_BM Wiggler | `Wiggler(E0 = 25, dE = 24, phase = 0, randomphase = 1, Ee = 2.4, Ie = 0.4, B = 1.6, K=10, Nper=41)`
MaxIV / Bloch | `Undulator(E0=0.6,dE=0.4,Ee=1.5,dEe=((6e-9)*(60e-12))/1.5,Ie=0.5,tbunch=43,K=5.6,gap=14e-3,Nper=187,lu=84e-3,sigey=1.3e-5,sigex=185e-5,sigepx=32e-6,sigepy=4.6e-6,focus_xw=1.1e-3,focus_yh=1.1e-3,dist=zm_mirror1,E1st=1.0018*E0/5)`
MaxIV DanMAX | `Undulator(E0=35, dE=0.05, E1st=E0/15, dist=20, Ie=0.5, Ee=3.0, dEe=0.0008, K=0, B=0,quick_integ=1, Nper=187, lu=0.016, sigex=53.66e-6, sigey=4.008e-6, sigepx=5.963e-6, sigepy=2.004e-6)`
SOLEIL U18 ANATOMIX (long) | `Undulator( E0=17, dE=1, Ee=2.75, dEe=0.001, Ie=0.5, K=1.03118, Nper=140, lu=32e-3, sigey=6.17e-6, sigex=0.29979e-3, sigepx=0.01226e-3, sigepy=1.1e-6, dist=50, E1st=12.400)`
SOLEIL U24 PX2a (medium straight section) | `Undulator( E0=12.65, dE=1, Ee=2.75, dEe=0.001, Ie=0.5, K=1.788, Nper=80, lu=24e-3, sigey=9.3e-6, sigex=215.7e-6, sigepx=29.3e-6, sigepy=4.2e-6, dist=29.5, E1st=12.400)`

In order to model the SOLEIL photon sources, we first need to refer to the storage ring parameters for SOLEIL:

- [https://www.synchrotron-soleil.fr/en/research/sources-and-accelerators/electron-beam-parameters/transverse-size-electron-beam-source](https://www.synchrotron-soleil.fr/en/research/sources-and-accelerators/electron-beam-parameters/transverse-size-electron-beam-source)

Then we may for instance look at the [LUCIA](https://www.synchrotron-soleil.fr/fr/lignes-de-lumiere/lucia) beam-line (SD03C) which is illuminated with an Undulator HU52 "Apple II" type (NdFeB magnets), 32 periods, gap 15-150mm, variable linear polarization, left and right circular polarizations, operating on harmonics 3 to 21. The energy range is 0.6-8 keV on LUCIA, and 0.35-2.5 keV on DEIMOS. 


HU52 parameter | symbol/unit | value
------------------------|-----------|-----
Period | (mm) | 52.4
Nb of periods  | lambda_U | 32
Gap | mm | 15.5-150
Field remanence | Br(T) | 1.26
Magnetic field Z, max | Bz (T) | 1.974*exp(-3.1754*Gap/lambda_U)
Magnetic field X, max | Bx (T) | 1.901*exp(-4.3387*Gap/lambda_U)
Function beta horizontal |	beta_x (m)	|1.4
Function beta vertical|	beta_z ( m)	|1.4
emittance horizontale | ex (pm.rad)	|82.0
betatron	coupling | % | 30.0
emittance vertical |	ez (pm.rad)	|24.6
dimension RMS horizontal|	sigma_x (µm)	|188(current)/10.7(upgrade)
dimension RMS vertical	|sigma_z (µm)	|8.2(current)/5.9	
divergence RMS horizontal	|sigma_x' (µrad)	|25.2(current)/7.7(upgrade)
divergence RMS vertical	|sigma_z' (µrad)|	6(current)/4.2(upgrade)

The corresponding HU52 Undulator component parameters are then:
``` c
Undulator(
  E0     = 3,
  dE     = 2.9,
  Ee     = 2.75,
  dEe    = 0.001,
  Ie     = 0.5,
  B      = 0.42, // for a 15.5 mm gap
  Nper   = 32,
  lu     = 52.4e-3,
  sigex  = 188e-6,
  sigey  = 8.2e-6,
  sigepx = 25.5e-6,
  sigepy = 6e-6) 
```

1. Start a new simulation and insert an Undulator source in it. The Undulator component has many possible parameters. In McXtrace, the (0,0,0)-point is taken to be the exit plane of the undulator. X is left-wise, Y is vertical, Z is forward. Use the typical HU52 Undulator component parameters for LUCIA or DEIMOS undulators.

2. Insert monitors **20 m** downstream: for instance one PSD "x y", one energy-resolved monitor, one divergence monitor "dx dy". Make sure that the monitors are big enough to catch all the radiation you expect, including the energy range. By using the `Monitor_nD`, you may add the "all auto" option to automatically adapt the monitor bounds to catch all photons. The corresponding code for the PSD could be `Monitor_nD(xwidth=0.1, yheight=0.1, bins=512, options="x y, all auto")` or the `PSD_monitor` component. You should get results such as ![HU52](images/mcplot_1.png?raw=true "")

3. Repeat the simulation with the expected SOLEIL-II Upgrade storage ring. Compare results in photon beam size and divergence.
 

References:

- F. Briquez et al., Proceedings of FEL08, Gyeongju, Korea, [https://accelconf.web.cern.ch/fel2008/papers/tupph015.pdf](https://accelconf.web.cern.ch/fel2008/papers/tupph015.pdf)
- T. Moreno et al., J Sync Rad 19 (2012) 179, [https://journals.iucr.org/s/issues/2012/02/00/kt5033/index.html](https://journals.iucr.org/s/issues/2012/02/00/kt5033/index.html)
- T. Moreno et al., [Proceedings Volume 8141, Advances in Computational Methods for X-Ray Optics II; 81410H (2011) DOI: 10.1117/12.893778](https://www.researchgate.net/publication/258548494_Undulator_emission_analysis_Comparison_between_measurements_and_simulations)
- M.E. Couperie 2013, [https://accelconf.web.cern.ch/ipac2013/talks/mozb102_talk.pdf](https://accelconf.web.cern.ch/ipac2013/talks/mozb102_talk.pdf)

## Exercise: Using the native McXtrace Bender model for SOLEIL photon sources

Let's now model a Bender at SOLEIL. For this we use the `Bending_magnet` component. Looking at its documentation, you will find that a typical use at SOLEIL is (for the ROCK bender):

``` c
Bending_magnet(
   E0 = 20, dE = 19, Ee = 2.75,
   Ie = 0.5, B = 1.72, sigey=9.3e-6, sigex=215.7e-6)
```
where we refer to the 'medium straight section' e-beam cross-section.

1. Create a copy of the above model, and change the Undulator for a Bending_magnet feeding the ROCK beam-line (E0=4.5-40 keV).

2. Update the e-beam parameters for the SOLEIL-II storage ring. Compare results in photon beam size and divergence.




---

## Exercise: Use MCPL files to couple to e.g. SRW and other codes (optional)

### Using the Undulator model

The typical Undulator parameters for a SPring-8 insertion device are:
``` c
    E0=13, dE=1, Ee=8, dEe=0.001, Ie=0.1, 
    K=1.03118, Nper=140, lu=3.2e-2, 
    sigey=6.17e-6, sigex=0.29979e-3, sigepx=0.01226e-3, sigepy=1.1e-6, 
    dist=20, E1st=12.400
```

1. Copy the initial SOLEIL Undulator model above, and use the SPring-8 parameters.

2. Perform a simulation using a large energy range, and look at the beam energy and spatial distribution, 20 m down-stream.

### Using MCPL files which store photon events

We will now use a different utility to drive a McXtrace-simulation: MCPL. 

[MCPL](https://mctools.github.io/mcpl/) is an interchange file format to communicate with e.g. GEANT4, PHITS, MCNP, and SRW. The SOLEIL Optics group is developing the [OptiX](https://gitlab.synchrotron-soleil.fr/OPTIQUE/optical-simulation/pyoptix) code which generates MCPL as well.

In this case the MCPL-file is actually generated using [SRW](https://www.github.com/ochubar/SRW). The expert tool for this purpose, a C++-program [`srw2mcpl`](https://github.com/McStasMcXtrace/srw2mcpl) is in development status and based on SRWlib and MCPL. The tool makes repeated calls to SRW and generates rays from that, and is not in its current form in an end-user state. 

The method of calling SRW pr. ray is rather slow, so for this tutorial (to save time) we provide a pre-generated MCPL-file that you may use. We are working on solutions to make the program available to McXtrace users in an easy-to use form. 

The file you need is called [sp8stdU.mcpl.gz](data/sp8stdU.mcpl.gz?raw=true ""). There is also a bigger version of this same file for better sampling [sp8stdUl.mcpl.gz](data/sp8stdUl.mcpl.gz?raw=true ""). But this obviously takes longer to download.

1. Use the McXtrace component **MCPL_input** to read rays from it and start them in a McXtrace simulation. 
2. Insert a downstream `PSD_monitor` and an `E_monitor` to catch the generated radiation. Leave some room (2 m or so) between the MCPL-file and the monitor. SRW considers the undulator centre its reference point and so rays may actually originate there. 
3. What was the fundamental energy of the 1st harmonic?
4. This procedure relies on the undulator spectrum being sufficiently sampled by the `srw2mcpl`-program. Determine the sampling limits of the file using your monitors. 
5. Compare with the pure McXtrace Undulator model above.


