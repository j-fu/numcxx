numcxx - a collection of python compatible classes for linear algebra in C++
======================================================================
Authors: 
    - Jürgen Fuhrmann (http://www.wias-berlin.de/~fuhrmann)
    - Timo Streckenbach (http://www.wias-berlin.de/~strecken)

## Rationale

- Provide  lightweight  multidimensional  linear algebra  classes  for
  C++11

- Keep the  code reasonably  simple and  transparent for  beginners in
  order to be useful for teaching

- Study expressive capabilities of modern C++ in the context of numerical
  algorithms

- Zero copy,  reference  counted exchange  of  array  data with  other
  packages without copying data, in particular
  - [LAPACK](http://www.netlib.org/lapack/) (to be installed via system installer)
  - [UMFPACK](http://faculty.cse.tamu.edu/davis/suitesparse.html) (to be installed via system installer)
  - [triangle](https://www.cs.cmu.edu/~quake/triangle.html) (provided with numcxx)
  - [iterative method templates](http://www.netlib.org/templates/cpp/) (provided with numcxx)

- interface to python/numpy (work in progress)


Many  ideas  used behind  this  library  have  been developed  in  the
framework of  the WIAS  [pdelib](http://pdelib.org) project,  and they
evolved over many  years.  C++11 allows for  a standard-conforming and
concise implementation, and so this code is as well a concept study.

Anyone  expecting a more  sophisticated package  may have  a look  at these
projects:

- [Eigen](http://eigen.tuxfamily.org)
- [Armadillo](http://arma.sourceforge.net/)
- [Blaze](https://bitbucket.org/blaze-lib/blaze/overview)
- [Trilinos/DOMI](https://trilinos.org/packages/domi)

## Further information
 - [Public mercurial repository (on bitbucket)](https://bitbucket.org/j-fu/numcxx)
 - [Doxygen documentation (via WIAS)](http://www.wias-berlin.de/people/fuhrmann/numcxx/html/index.html)
 - Introduction: [(local)](intro.md) [(via WIAS)](http://www.wias-berlin.de/people/fuhrmann/numcxx/html/md_doc_intro.html) 
 - Installation: [(local)](installation.md) [(via WIAS)](http://www.wias-berlin.de/people/fuhrmann/numcxx/html/md_doc_installation.html) 

