# Lecture on some advanced tricks you can pull since we have low level access to the code.

Also includes some advanced grammar features like **EXTEND**, **GROUP**, and **JUMP**
1. Focal point monitor using `Arm+Monitor_nD`
2. Model scintillator efficiency using `Filter+PSD`
3. Write your own component from scratch e.g. the semi-permeable mirror.
4. Union?


# COMPONENT syntax in beam-line models

```
TRACE
...

{SPLIT N} COMPONENT name = comp(parameters) 
  {WHEN condition}
  AT (...) [RELATIVE [reference|PREVIOUS] | ABSOLUTE]
  {ROTATED {RELATIVE [reference|PREVIOUS] | ABSOLUTE} }
  {GROUP group_name}
  {EXTEND %{ C_code %} }
  {JUMP [reference|PREVIOUS|MYSELF|NEXT] [ITERATE number_of_times
| WHEN (condition)] }

...
END
```
