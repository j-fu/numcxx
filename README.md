numcxx - a collection of python compatible classes for linear algebra in C++ {#mainpage}
======================================================================
Authors: 
    - Jürgen Fuhrmann (http://www.wias-berlin.de/~fuhrmann)
    - Timo Streckenbach (http://www.wias-berlin.de/~strecken)

## Rationale

- Provide  lightweight  multidimensional  linear algebra  classes  for
  C++11

- Keep the  code reasonably  simple and  transparent for  beginners in
  order to be useful for teaching

- Study expressive capabilities of modern C++

- Efficient,  reference  counted exchange  of  array  data with  other
  packages without copying data
  - LAPACK
  - UMFPACK
  - triangle (soon)

- SWIG based interface to python/numpy


Many  ideas  used behind  this  library  have  been developed  in  the
framework of  the WIAS  [pdelib](http://pdelib.org) project,  and they
evolved over many  years.  C++11 allows for  a standard-conforming and
concise implementation, and so this code is as well a concept study.

Anyone  expecting more  sophisticated work  may have  a look  at these
projects:

- eigen
- armadillo
- Trilinos/DOMI
- marray
- numcpp

