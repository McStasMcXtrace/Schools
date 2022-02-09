Build along session where we will build an instrument from scratch

In this exercise we will all together build an extremely simple beamline consisting only of a source, two slits, and two different detectors/monitors.
# Step 1. Open mxtrace.
Start by either typing mxgui in ta terminal, double clicking the McXtrace icon, or choosing McXtrace from the application/start menu. Whatever is approprioate for your system,
The McXtrace main window should now open, looking like this.
![McXtrace main window](images/mxgui_main.png?raw=true "")

# Step 2. Start a new instrument description.
In the 'File' menu, please choose the 'New Instrument' item. This opens a file dialog window where you can choose a name for your new instrument file. Do so and click save. This now opens the Editor window of the GUI where you can edit a very simple template instrument file. 

![McXtrace editor window](images/mxgui_editor.png?raw=true "")

A McXtrace instrument description is a text-file just like any other program you write. Generally the suffix `.instr` is used but you could use whatever you like for this.

# Step 3. Insert a source.
The editor can help you write instrument files, but really is just an editor. If you prefer to use something else this works equally fine.
Some of the allowed McXtrace code blocks have already been defined in the templat file - Position your cursor below the keyowrd **TRACE** and go to the Insert menu. Choose the item `Sources -> Source_pt`. This opens the insert component dialog.

# Step 4. Fill in some source parameters.
Enter some values into the boxes on the right hand side of the window. For instance you may choose: 
- focus_xw=1e-3
- focus_yh=1e-4
- dist=2
- E0=15
- dE=1

Lastly click 'Insert'. Ths inserts code at the cursor position corresponding to the parameters which you entered.
We could just as easily have entered these details by hand, in which case the editor helps us a little bit by autocompleting McXtrace keywords.

The Insert dialog is not clever enough to figure out where your code should go, so you _do_ have to position the cursor where the component should be inserted.

# Step 5. Insert a Monitor.
Now that we have generated some photons, we shall also try to detect them. Position the cursor below the source component from before and go to `Insert -> Monitor -> PSD_monitor`. This is a very simple model of a flat area detector divided into nx by ny pixels. Please enter dimensions of the monitor you want and also a filename in which it should save results.  Further, we need to place the monitor downstream of the source. Do this by entering a value (e.g. 2) in the 3rd (for Z coordinate) box of the AT row in the bottom of the insert dailog window. Remember that in McXtrace, Z generally points towards the next component. Most of the time this coincides with the optical axis.
We now have a complete (albeit boring) beamline which we can run. Click the `Run`-button in the main GUI to do so.

You should now be presented with a new dialog asking you to set some parameters for your simulation. For now the default should suffice, so click the `Run` button to start the simulation. If all  goes well, this should return us to the main window where we may click `Plot` once the simulation is finished.
This ought to give you a plot like:

![First Simulation Results](images/mxplot_first_sim.png?raw=true "")

# Step 6. Get acquainted with the documentation tools
From the mxgui Help menu, the following entry points are available

![mxgui Help  Menu](images/mxdocfromGUI.png?raw=true "")

mxdoc current instrument (will present you with information on the currently loaded instrument, if available)

mxdoc Component Reference (will give you access to online-docs for all installed McXtrace components)

![mxdoc Component Reference](images/mxdoc_browser_overview.png?raw=true "")

![mxdoc component view](images/mxdoc_component.png?raw=true "")

mcxtrace User Manual (will open your McXtrace user manual as a PDF)

![mxdoc User Manual](images/mxdoc_manual.png?raw=true "")

mcxtrace Component Manual (will open your McXtrace component manual as a PDF)

![mxdoc User Manual](images/component_manual_front.png?raw=true "")

mcxtrace Web Page (will take you to the McXtrace website)

