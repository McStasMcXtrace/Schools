# Solutions and ideas
The instruments in this "solution" folder are:
* [`mono.instr`](mono.instr) is the result of the exercise steps
* [`mono_optims.instr`](mono_optims.instr) is meant as an inspiration catalogue for a subset of many optimistaions possible. 

The listing below highlights the differences between the two instruments:

```diff
2c2
< * Instrument: mono
---
> * Instrument: mono_optimscd
13c13,15
< * Instrument "solution" for Monochromator exercise
---
> * Instrument "solution" for Monochromator exercise, with further optimisation:
> * 1) A parameter lambda - if set will calculate A1 and A2 for this wavelength
> *    (Otherwise lambda is calculated from the A1 setting)
14a17,26
> * 2) A safety check that lambda is in source_lam_min/max interval
> * 
> * 3) Calculation of "focusing / roi" to not overshoot mono cross-section too much
> *
> * 4) Mono_flat -> Mono_curved with vertical focusing option
> *
> * 5) Source-mono collimation can be controlled by soller collimator if coll_in is set
> *
> * 6) Mono-sample collimation can be controlled by soller collimator if coll_out is set
> *
21a34,37
> * lambda:         [AA] Automatic setup of A1/A2 for this wl if set
> * coll_in:    [arcmin] If set, defines in-pile horizontal collimation by soller
> * coll_out:   [arcmin] If set, defines post mono horizontal collimation by soller
> * RV:            [deg] Vertial focusing of monochromator
28c44
< DEFINE INSTRUMENT mono(source_lam_min=0.5, source_lam_max=6.5, A1=17.3472,A2=34.6943)
---
> DEFINE INSTRUMENT mono_optims(source_lam_min=1.95, source_lam_max=2.05, A1=17.3472,A2=34.6943, lambda=2, coll_in=30, coll_out=30, RV=0.6)
31a48,52
>   double Theta; // Internal variable for Mono bragg condition calc
>   double Mono_D; // Mono d-spacing
>   double Mono_w; // Mono width
>   double Mono_w_src; // Mono @A1, width seen from source
>   double MSdist; // Mono-sample distance
39a61,87
>   // Initalise constant(s)
>   Mono_D=3.3539; // Mono PG 002
>   Mono_w=0.05;    // Mono width
>   MSdist=2.1;    // Mono-sample distance
> 
>   /* Optimisation 1) - 
>      Check if we are provided with a "target" wavelength */
>   if (lambda) {
>     Theta = RAD2DEG * asin(lambda/(2*Mono_D));
>     A1 = Theta;
>     A2 = 2*Theta;
>     printf("Input (target) lambda defined, setting A1=%g, A2=%g\n",A1,A2);
>   } else {
>     Theta = A1;
>     lambda = 2*Mono_D*sin(DEG2RAD*Theta);
>     printf("Target lambda calculated from A1 setting, lambda=%g\n",lambda);
>   }
> 
>   /* Optimisation 2) */
>   if(lambda < source_lam_min || lambda > source_lam_max) {
>     fprintf(stderr,"ERROR: Mono wavelength %g is not in interval from source [%g %g]\n",lambda,source_lam_min,source_lam_max);
>     exit(-1);
>   }
> 
>   /* Optimisation 3) - add 5 mm overshoot */
>   Mono_w_src = 0.005 + Mono_w*sin(RAD2DEG*A1);
> 
60a109,118
> /* Optimisation 5) */
> COMPONENT InpileColl = Collimator_linear(
>   xmin = -0.05,
>   xmax = 0.05,
>   ymin = -0.05,
>   ymax = 0.05,
>   length = 0.10,
>   divergence = coll_in) 
> WHEN (coll_in) AT (0, 0, 3.5) RELATIVE  Source
> 
65c123,125
< COMPONENT monochromator_flat = Monochromator_flat()
---
> /* SPLIT added to boost our stats beyond mono 10-fold */
> /* + optimisation 4) */
> SPLIT COMPONENT monochromator_flat = Monochromator_curved(zwidth=Mono_w,RV=RV,NH=1)
71a132,141
> /* Optimisation 6) */
> COMPONENT MonoColl = Collimator_linear(
>   xmin = -0.05,
>   xmax = 0.05,
>   ymin = -0.05,
>   ymax = 0.05,
>   length = 0.10,
>   divergence = coll_out) 
> WHEN (coll_out) AT (0, 0, MSdist/2.0) RELATIVE  Mono_out
>   
76,77c146,147
<     xwidth=0.05,yheight=0.05)
< AT (0, 0, 2) RELATIVE Mono_out
---
>     xwidth=0.1,yheight=0.1)
> AT (0, 0, MSdist-0.1) RELATIVE Mono_out
79c149
< COMPONENT PSD = PSD_monitor(xwidth=0.05,yheight=0.05)
---
> COMPONENT PSD = PSD_monitor(xwidth=0.1,yheight=0.1)
83c153
< AT (0,0,2.1) RELATIVE Mono_out
---
> AT (0,0,MSdist) RELATIVE Mono_out
```
