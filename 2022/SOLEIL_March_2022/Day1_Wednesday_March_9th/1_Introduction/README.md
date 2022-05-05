# McXtrace training: Introduction

This is an introductory lecture on  [McXtrace](http://www.mcxtrace.org). We present how Monte-Carlo ray-tracing works, how  [McXtrace](http://www.mcxtrace.org) has been designed, the notions of *components* and *instruments*, and actually how to get help.

## Introduction to Monte-Carlo techniques

The  [McXtrace](http://www.mcxtrace.org) software is a ray-tracing Monte-Carlo code, that is it traces rays (photons) and involves random numbers in their generation and/or interaction with matter and beam-line parts.

The Monte-Carlo technique consists in sampling randomly a phase-space (the space of variables that describe a system), and compute integrated quantities (macroscopic observable).

If we sample a variable *x* with *N* random values in a range `[a,b] `, then the mean distance between two adjacent such values converges to (for large *N*):

<dl><center>
&delta;x → (b-a)/N 
</center>
</dl>

and then the integral can be approximated (for large *N*):

<dl><center>
∫<sub>[a,b]</sub> f(x) dx → (b-a)/N  ∑ f(x) 
</center>
</dl>

Surprisingly, it is possible to demonstrates that the Monte-Carlo integration approximation converges faster than any other method (trapezoidal, Simpson, etc) when the dimensionality *d* of the variable space is larger than 5.

For instance, with particles such as photons and neutrons, the typical dimensionality is 10, which ensures accurate integral estimates with few random events. More specifically, 10<sup>3</sup> and  10<sup>5</sup> random events in the integral provide an accuracy of about 10% and 1% respectively.

### Improving the efficiency of the Monte-Carlo methods

In practice, there are well known techniques to speed-up the random event integral computations. In many cases, the random events are attached a probability *weight*. For instance with photons, this can be considered as a sort of particle pack (kind of a small flux fraction). 

When it is possible to evaluate numerically the probability that an event occurs, we can just multiply the weight by this probability and progress to the next event. This avoids casting a random number and becomes a deterministic implementation.

More generally, when the probability of an event follows a non flat distribution (non equiprobable choices) there are values of the random variable that contribute 'more' to the integral. Then, it makes sense to count more of these, and less in the low-probability areas. We can then consider this as a *variable change* from a flat to a non-flat distribution. Our tuned random distribution can be used as is, and we need to correct for the associated Jacobian. This widely used technique is called *importance sampling*. 

An other technique consists in identifying a sub-region of the phase-space that contains more significant contributions to the statistics, and in which further processes involve random numbers. Then, we can treat the inner and outer regions separately, as long as there are enough events in the inner sub-region and the number of random casts that follow is of the order of the system dimensionality *d*. This is called *stratified sampling'*.

References:
1. F. James. In: *Rep. Prog. Phys.* **43** (1980), p. 1145.

## History, construction of McXtrace

The  [McXtrace](http://www.mcxtrace.org) software has been derived in 2009 from the sister project [McStas](http://www.mcstas.org) for neutron ray-trace (which dates from 1998). It shares the same core infrastructure: the way to write beam-line descriptions, the way to write components (parts used in beam-lines), the "grammar" 

## "Components" and "Instruments"

## Tips and Tricks / How to get help

---
*McXtrace training - 2022*