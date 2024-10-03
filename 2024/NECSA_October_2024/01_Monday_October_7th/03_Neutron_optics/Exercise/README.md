## Guides
The exercise starts with the [`Exercise_guides.instr`](Exercise_guides.instr) file, it contains a source, straight guide and relevant monitors. It already has an instrument parameter for controlling the length of the guide.
Solutions are available in `solution/Exercise_solution_straight.instr` and `solution/Exercise_solution_elliptic.instr`.

### TASKS
- Task 1
    - Compare output for two different guide lengths (e.g. triple the length)
- Task 2
    - Introduce a parameter that control width of the guide (default 3 cm)
    - Compare two runs with different guide widths (e.g. half and double the width)
- Task 3
    - Check how much gravity impacts the output (set `G` to 0 and Jupyter `2.528*G`)
- Task 4 (optional)
    - Exchange the last 20% of the guide with an elliptic nose.
    - See the geometry with `mcdisplay` (3D view)
    - Identify how the resulting beam has changed

### HINTS
* Use a scan in McStas to run two simulations with different guide lengths and compare
    * Set steps to 2 and in the field for `guide_length` have for example 5,15
    * The display will show the evolution in intensity on the monitor
* When adding a guide_width, remember to adjust the focusing of the source
* To add the elliptic nose use Elliptic_guide_gravity
* ```mcdoc Elliptic_guide_gravity```

### INTERPRETATION
The geometry of the guide heavily impacts the quality of the delivered beam.

### EXTRA
Investigate how the m value impacts the beam delivered by the guide.


