# Lecture on samples, plus includes advanced grammar for the instrument

Also includes some advanced grammar features like **EXTEND**, **GROUP**, and **JUMP**
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
