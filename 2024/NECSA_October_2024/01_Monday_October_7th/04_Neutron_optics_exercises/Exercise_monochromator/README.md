## Monochromator

First a briefing by Peter on McStas monochromator models:
[monochromator slides](Monochromators.pdf)

## TASKS

(If you get stuck along the way, there is a solution [here](solution/mono.instr))

1. Create yourself a new instrument file and name to your liking
2. Add a `Source_gen` to your instrument, by copy-paste from:

```
COMPONENT  Source = Source_gen(
    radius = 0.0905, dist = 4, focus_xw = 0.1, focus_yh = 0.1,
    Lmin = source_lam_min, Lmax = source_lam_max, I1 = 0)
  AT (0, 0, 0) RELATIVE origin
```
3. Add corresponding instrument input paramters for your source: `source_lam_min=0.5, source_lam_max=6.5`,
4. Insert an `L_monitor` a short distance after the `Source` to measure the initial spectrum, look up reasonable parameters via the component insertion dialogue.
5. Run a simulation to verify output on the. monitor
6. Next we need to add two `Arm` components for the positioning of our Monochromator:
   * Instance `Mono_arm` with `AT (0,0,4) RELATIVE Source` and `ROTATED (0,A1,0) RELATIVE Source`
   * Instance `Mono_out` with `AT (0,0,4) RELATIVE Source` and `ROTATED (0,A2,0) RELATIVE Source`
7. Add input parameters `A1=0, A2=0` for control of the geometry of your instrument
8. Use the `mcdoc` utility (press `Docs` on `mcgui`) to read up on the component `Monochromator_flat`
9. Now insert a Monochromator_flat between `Mono_arm` and `Mono_out` using (all default) reflection of PG
10. Insert another `L_monitor` instance 2 m after the `Mono_out` instance to (once we are ready) measure the reflected beam from the monochromator. Use the same range of wavelength measurement as the earlier wavelength monitor.
11. Following 0.001 mm after the (10.) L_monitor, insert a 0.1 x 0.1 m `PSD` to look at the reflected beam spot.
12. Verify that you have built something reasonable by means of a `TRACE` setting `A1=45` and `A2=90`, it should all togetherlook something like this:
![3D](pics/mcdisplay_1.png)
![3D](pics/mcdisplay_2.png)
