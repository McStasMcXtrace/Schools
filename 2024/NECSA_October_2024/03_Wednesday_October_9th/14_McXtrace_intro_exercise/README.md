# McXtrace introductory exercise.

![McXtrace](../../pics/mcxtrace-logo.png  "McXtrace")

In this exercise we shall build an extremely simple beam-line consisting only of a source, a slit, and a detector/monitor.

## Step 1. Open McXtrace.

Start by either typing `mxgui` in a terminal, double clicking the McXtrace icon, or choosing McXtrace from the application/start menu.
The McXtrace main window should now open, looking like this.

![McXtrace main window](pics/mxgui_main.png?raw=true "")

## Step 2. Start a new instrument description.

In the 'File' menu, please choose the 'New Instrument' item. This opens a file dialog window where you can choose a name for your new instrument file (e.g. `test.instr`). Do so and click _Save_. 

A McXtrace instrument description is a text-file just like any other program you write. Generally the suffix `.instr` is used but you could use whatever you like for this.

This opens the Editor window of the GUI where you can edit a very simple template instrument file. You may equally open the model with any text editor.

:warning: There is a **bug** with the 3.5.1 `mxgui`. The initial opening of the new file may break and exit `mxgui`. Restart `mxgui` and reopen the new file. An external editor window will pop-up (e.g. `gedit`). The ability to insert pre-defined components is not available.

![McXtrace editor window](pics/mxgui_editor.png?raw=true "")

## Step 3. Insert a source.

The editor can help you write instrument files, but really is just an editor. 
If you prefer to use something else, this works equally fine.

Some of the allowed McXtrace code blocks have already been defined in the 
template file. You can identify the usual sections and keywords:

- a header that describes the purpose of the model. Please keep information up-to date, especially the input parameter details.
- the `DEFINE INSTRUMENT` line that gives the name of the model, and its input parameters.
- the `DECLARE` section where C-style variables may be listed, and may be used elsewhere.
- the `USERVARS` section where to define C-style variables attached to each particle.
- the `INITIALIZE` section which contains C-code to execute before shooting the first particle.
- the `TRACE` section which contains a ordered list of `COMPONENT` lines.
- the `FINALLY` section which contains C-code to execute at the end.
- the `END` keyword.

Only the  `DEFINE INSTRUMENT`, `TRACE` and `END` parts are mandatory.

Position your cursor below the keyword **TRACE** and insert the code:
```c
COMPONENT src = Source_pt()
AT (0,0,0) ABSOLUTE
```

This way we have positioned a point source called `src` (this name is fully yours), for which some parameters must be defined. Not all parameters are required, and often there are different ways to specify the same information.
Component parameters are usually a set of:

- geometrical (width, height, shape, ...)
- physical (energy, density, ...)
- computational (probabilities, statistical focusing)
- input/output (file names, ...)

It is recommended to open the documentation by clicking the 'Docs' button in the main interface. A browser window will show-up.

![mxdoc Component Reference](pics/mxdoc_browser_overview.png?raw=true "")


## Step 4. Fill in some source parameters.

It is first required to get documentation about the components to be inserted. 

![mxdoc component view](pics/mxdoc_component.png?raw=true "")

Most components contain an `Example:` line which provides a basic syntax with main parameters.

Select the component `Source_pt` from the component list (browser). This is a point source, which illumination is isotropic.
Looking at the `Source_pt` documentation, we select:

- focus_xw=1e-4
- focus_yh=1e-4
- dist=2
- E0=15
- dE=1

so that the final component description should be:

```c
COMPONENT src = Source_pt(focus_xw=1e-4, focus_yh=1e-4, dist=2, E0=15, dE=1)
AT (0,0,0) ABSOLUTE
```

These settings mean that the beam size will be 1x1 mm<sup>2</sup> at 2 m distance (along `Z`), and the beam energy will be centred around 15 keV.

## Step 5. Insert a Slit and a Monitor.

Now that we have generated some photons, we shall shape the beam, and measure the spatial distribution. 

With the same methodology, insert a `Slit` component 0.7x0.7 mm<sup>2</sup> at 2 m, and then  a `PSD_monitor` component at 5 m. 
Set the monitor width and height so that the beam spot spreads on most of the surface.

We now have a complete (albeit too simple) model which we can run. 
Save the file and click on the `Run`-button in the main GUI.

You should now be presented with a new dialog asking you to set some parameters for your simulation. For now the default should suffice, so click the `Run` button to start the simulation. If all goes well, this should return us to the main window where we may click `Plot` once the simulation is finished.
This ought to give you a plot like:

![First Simulation Results](pics/mxplot_first_sim.png?raw=true "")

## Step 6. Add more monitors and parameterise the model.

It is possible to change the monitor type (and its parameters), or add others.

We suggest that you add an `DivE_monitor`, i.e. a monitor which first axis is 
Energy and 2nd is the horizontal divergence. Position this new component a few cm 
away after the previous one, and set some reasonable parameters.

**Question:** how would you get an energy/vertical divergence monitor ? (tip: `ROTATED`).

Now, we wish to be able to change the source energy without changing the model 
description file. The `DEFINE INSTRUMENT template_simple(parameter1=1)` is there 
for that purpose.
Change the default `parameter1=1` (which is not used here), by e.g.:
```c
E0=15
```
This syntax sets an `E0` parameter with a default value of 15.
You may at the same time change the name of the model to match your file name 
(e.g. `test.instr` -> "DEFINE INSTRUMENT **test**(".

Now use the `E0` model parameter in the `Source_pt` parameters. It could even be 
used for the `DivE_monitor(Emin=..., Emax=...)`, so that the monitor adapts its 
range as a function of the source.

Save and run the model again, and look at the results by clicking the Plot button.
As you can now see, the input model parameter is now `E0`, set to 15 at start, 
but allowing to change its value without touching the model description.

The two generated files from the monitors are text based.
