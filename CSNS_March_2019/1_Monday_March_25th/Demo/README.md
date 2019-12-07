## A quick McStas demonstration

### Creating a new instrument

1. Open mcgui on your computer
2. Choose ```File ->  New Instrument``` and save under a relevant location / filename
![Create new instrument from template](images/new_instrument.png?raw=true "")

### Insert a source

3. Position the cursor below the comment ```// insert components here (e.g. Insert -> Source -> ...)``` 
![Insert the source](images/insert_source.png?raw=true "")
4. Use ```Insert->Source->Source_simple``` and fill in
   * a radius of 0.1 m (also default)
   * a dist of 10 m, ``` focus_xw x focus_yh ``` of 0.01 x 0.07 m
   * choose a ```lambda0 ``` of 2 Å and ``` dlambda```  of 0.02 Å
![Choose the source parameters](images/source_parameters.png?raw=true "")
5. Press ``` Insert ``` and check that the input looks good
![Inspect inserted code](images/source_code_now_in_editor.png?raw=true "")

### Insert a sample

6. Similarly, use ```Insert->Sample->Powder_N``` to insert a sample after the source
![Insert the sample](images/insert_sample.png?raw=true "")
7. Use ```"YBaCuO.lau"``` as reflections, radius of 0.005 m and yheight of 0.07 m
8. Place at 10 m from the source
![Choose the sample parameters](images/sample_parameters.png?raw=true "")
![Inspect inserted code](images/sample_code_now_in_editor.png?raw=true "")

### Insert a PSD

9. Use ```Insert->Monitor->PSD_monitor``` to insert a position-sensitive monitor
![Insert a PSD](images/insert_PSD.png?raw=true "")
10.  Fill in ```"PSD.dat"``` as filename and chose ```width``` and  ```width``` both at 2 m, ```nx``` and ```ny``` both set to 401
11. Place the PSD 1m after the sample
![Choose the PSD parameters](images/PSD_parameters.png?raw=true "")
![Inspect inserted code](images/PSD_code_now_in_editor.png?raw=true "")

### Run your instrument - trace mode

12. Press run on the main interface, if all is good you should now see this dialogue, where you should chose "Trace" at the top
![Run in Trace mode](images/run_mcdisplay.png?raw=true "")
13. Press ```Start``` and a browser should appear
![Trace in browser](images/mcdisplay_output.png?raw=true "")
14. Press run again and choose Simulation, then Start
![Run in Simulation mode](images/run_simulation.png?raw=true "")
15. Pressing Plot will give you the data from the PSD
![Linear scale plot](images/lin_plot.png?raw=true "")
16. Pressing ```l``` (L) for log should give you visible powder lines
![Log scale plot](images/log_plot.png?raw=true "")
