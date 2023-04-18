# Exercise 2. Kirkpatrick Baez (KB) mirrors
The aim of this exercise is to replace CRL stack with a Kirkpatrick Baez focusing system in McXtrace. We will use two instances of a single curving mirror model: `Mirror_curved`, to focus
the beam in the same manner as exercise 1, i.e. focus the beam at a pt 10 m
downstream from the focusing system at 31.5 m downstream from the source.

<img src="https://www.j-tec.co.jp/english/wp-content/themes/j-tec-corp/images/page/img_en_kb-mirror3.jpg">

`Mirror_curved` is by default set in the YZ-plane. Think of how you can set your
mirror to deflect the beam in the desired angle. Remember that for a mirror
at a glancing angle &theta; degrees, the beam is redirected 2&theta;. The &theta; value may be an input parameter of the model, and its value is generally small (say 3 mrad).

## Reflectivity Data files

An extensive documentation on how to get and generate reflectivity data files for e.g. mirrors is accessible at:
- [HOWTO: McXtrace reflectivity files](https://github.com/McStasMcXtrace/McCode/wiki/HOWTO%3A-McXtrace-reflectivity-files) (requires IDL)
- [HOWTO: McXtrace reflectivity files XrayDB python library](https://github.com/McStasMcXtrace/McCode/wiki/HOWTO%3A-McXtrace-reflectivity-files-XrayDB-python-library)

These files are used by all mirrors.

## 1. Insert the 1st mirror
Start by re-using the source as defined in the [compound refractive lenses](../1_CRLs/README.md) exercise.

Insert an horizontal instance of a curved mirror at 31.5 mm. According to [KB 1948] the radius of curvature should be given by:

f_m = R sin(&theta;/2). 

Choose sensible length and width for your mirror and figure out what the glancing angle should be, and how curved it needs to be.
Input these parameters to the mirror component and set the AT to something meaningful. Also set a rotation to &theta; deg. around the second axis.

## 2. Insert an Arm to track the optical axis
In McXtrace, an `Arm`-component can be used to target an orientation change. An Arm does not do anything, it merely serves as a coordinate system, and so can be placed anywhere in your simulation.
Rotate the arm by &theta; RELATIVE to the mirror such that the z-axis again points along the optical axis.

This component helps to position further elements along the reflected beam.

## 3. Insert the 2nd mirror some distance behind the 1st
The 2nd mirror needs to be perpendicular to the first. Therefore we add an extra Arm, *a2*, in front of the mirror which is ROTATEd by 90 deg. around the Z-axis, then place the 2nd mirror RELATIVE to *a2*, and ROTATE it by &theta; around the Y-axis. The distance between the two KB mirrors is usually rather small, but you must avoid overlapping them.

Lastly, place another Arm after the 2nd mirror (but AT the same place) which is ROTATED another &theta; around the Y-axis just as we did for mirror 1. 

For the second mirror the curvature may be found (again following [KB 1948]) from

f_s = R/2 sin(&theta;).

## 4. Verify that the focal point is in the right place.
Just as before, use a monitor and scan it along the optical axis to find the focal point, _or_ scan a pinhole (Slit) along the axis and find the place where its throughput is maximal.

We suggest that you define an input parameter of the model (distance from the KB to the monitor/slit) so that you can scan it.


[KB 1948]: https://www.ncbi.nlm.nih.gov/pubmed/18883922
Kirkpatrick, P. and Baez, V. A., "Formation of optical images by X-rays.", J. O. Soc. A., 1948
