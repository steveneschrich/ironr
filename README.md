# ironr
A port of of IRON normalization to R.

This is a port of C code from libaffy into the R framework. Some of the principles
of this particular project include:

- **R style**. Where possible, C idioms will be replaced by R idioms so that the
 interface and usage will be as R-like as possible.
- **Equivalence vs. compatability**. The original software developed over a long period
of time. Where possible we will reframe certain implementations with equivalent
approaches.


## Major Functions

- iron.cel.findmedian() - most similar to libaffy findmedian program
- iron.cel.medoid() - analagous to findmedian but using more traditional metrics
