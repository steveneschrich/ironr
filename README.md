# ironr
A port of libaffy (primarily for IRON) normalization to R.

This is a port of C code from libaffy into the R framework. Some of the principles
of this particular project include:

- Minimizing dependencies. The original code had minimal dependencies for broad
 use, the goal would be (where possible) to do the same here.
- R style. Where possible, C idioms will be replaced by R idioms so that the
 interface and usage will be as R-like as possible.
 
Note: The goal is to use this reimplementation as the basis for a broader implementation
using Rcpp.
