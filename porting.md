# Porting IRON to R

This is a documenting outlining the thinking and steps associated with porting
the IRON algorithm into R as linked C object code. Currently, IRON lives within
a libaffy codebase (https://gitlab.moffitt.usf.edu:8000/WelshEA/libaffy).The
goal is to extract the C components that are computational into compiled R
library code (still in C) with necessary R wrappers to move data around. Integral
to this effort is the desire to use affyio, affy, oligo and other Bioconductor
packages as appropriate to reduce the amount of code as well as integrate better
with Bioconductor.

## libaffy
The libaffy library is a C code base for manipulating Affymetrix data. Eric Welsh
implemented his IRON algorithm within this code base. Therefore one of the first
steps should be to evaluate what code must be ported to R.

- libutils: Utility functions to make the library cross-platform compatible. Likely
  not necessary or implemented within R/C library of R.
- affyapps: These are the actual executables that call functionality within libaffy. As
  executable wrappers, these represent the API's for library functionality (merged
  with library parameters). 
  - findmedian
  - iron
  - iron_generic
  - pairgen
- libaffy
  - mas5/iron_bg.c  (not used, ignore)
  - mas5/iron_norm.c (the heart of IRON)
  - utils/ (some functionality as needed for IRON)
  - iron_generic/
  
