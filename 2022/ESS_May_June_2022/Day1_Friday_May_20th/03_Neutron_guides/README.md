## Guides
The exercise starts with the Exercise_guides.instr file, it contains a source, straight guide and relevant monitors. It already has an instrument parameter for controlling the length of the guide.
Solutions are available in Exercise_solution_straight.instr and Exercise_solution_elliptic.instr.

### TASKS
- Task 1
    - Compare output for two different guide lengths
- Task 2
    - Introduce a parameter that control width of the guide
    - Compare two runs with different guide widths
- Task 3
    - Check how much gravity impacts the output
- Task 4
    - Exchange the last 20% of the guide with an elliptic nose.
    - See the geometry with mcdisplay
    - Identify how the resulting beam have changed

### HINTS
* Use a scan in McStas to run two simulations with different input and compare
    * Set steps to 2 and in the field for guide_length have for example 5,15
    * The display will show the difference in intensity on the monitor
    * In display hold meta key and click on a graph to see the image comparison
* When adding a guide_width, remember to adjust the focusing of the source
* To add the elliptic nose use Elliptic_guide_gravity
* ```mcdoc Elliptic_guide_gravity```

### INTERPRETATION
The geometry of the guide heavily impacts the quality of the delivered beam.

### EXTRA
Investigate how the m value impacts the beam delivered by the guide.


