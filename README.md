numcxx - a collection of python compatible classes for linear algebra in C++ {#mainpage}
======================================================================

## Rationale
- Provide lightweight multidimensional linear algebra classes for C++11

- Keep the  code reasonably simple and transparent  for beginners in order  to be
  useful for teaching
     - Stay away from expression templates
     - Study expressive capabilities of modern C++

- Efficient, reference counted exchange of array data with other packages without copying
  data
     - python/numpy
     - triangle, TetGen etc.
     - LAPACK
     - UMFPACK

- SWIG based interface to python/numpy


Many ideas used behind this  library have been developed together with
or even  inspired by Timo  Streckenbach from the WIAS  pdelib project,
and   they   evolved   over   may   years.    C++11   allows   for   a
standard-conforming and concise implementation.

Anyone expecting more sophisticated work (expression templates in particular)
may have a look at these projects:

- eigen
- armadillo
- Trilinos/DOMI
- marray
- numcpp

