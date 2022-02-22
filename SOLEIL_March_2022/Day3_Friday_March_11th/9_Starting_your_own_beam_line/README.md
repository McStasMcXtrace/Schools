# Starting your own beam-line model: we help you !

In this 'anarchy' session, we shall do our best to get you started with your own beam-line model.

The methodology is rather simple:

1. Assemble a succinct list of elements that compose your beam-line. As a start, focus on the main components.
2. Identify which [McXtrace components](http://mcxtrace.org/download/components/) best match your beam-line elements. Prefer simple components when there are many possibilities.
3. Search for existing models that use these components. They are often in the 'Tests' category, but you may as well find a model which is pretty close in spirit to your beam-line. These examples will demonstrate how to actually use the components (syntax, parameters, etc) you envisage. Get inspiration here.
4. For each identified [McXtrace component](http://mcxtrace.org/download/components/), look at their main properties (parameters, geometry) and search for these in real life. This mean you will need to estimate plausible values, and potentially refer to some documentation. As a start, use approximate values.
5. Position these elements in space w.r.t. the other elements.
6. Once done, you may identify which values/parameters could be tunable in the beam-line model in order for instance to further scan or optimize their values. These can then be moved to the DEFINE INSTRUMENT line, and their variable be transferred where appropriate in the INITIALIZE and TRACE sections. set default values for these parameters.
7. Complete the header at the beginning of the .instr file to reflect your institution (SOLEIL), the beam-line, you as author, a description of the beam-line (copy-paste from [SOLEIL website](https://www.synchrotron-soleil.fr/en/beamlines)), a description of the model parameters (those from the step above).
7. Send us your beam-line, even basic/incomplete, so that we can further help you and enrich the gallery of SOLEIL beam-lines given as examples.

:bulb A special note on mirrors and monochromators: the positioning/orientation of these can be tricky. It is usually a good advice to search for examples that contain them, and copy/paste/adapt the corresponding code in your own model. You will save a lot of time scratching your head this way.

Of course, your model will hardly compile initially. Read the compilation output to decode the cryptic messages and identify where errors could be. An other strategy is to comment (use `/* ... */`) most of the components and leave only the source. Then compile and iteratively un-comment components one-by-one with intermediate compilation checks.

In the end, the structure of the beam line should have:

1. A photon source
2. Some optics (e.g. slits, monochromators, mirrors, lenses, ...)
3. A sample model as it provides a result which goes beyond a simple photon beam. The `PowderN` `Single-crystal` or `Absorption_sample` are good bets.
4. Monitors (you may position as many as required, not only at the end). The `Monitor_nD` is extremely versatile, but often you will start with e.g. `E_monitor` or `PSD_monitor`.


