numcxx - a collection of python compatible classes for numerics in C++
======================================================================

One of  the main goals: use  in teaching. So transparency  of the code
for people new to C++ is  important. Therefore numcxx does not provide
expression templates.   OTOH, based on  SWIG, numcxx provided  a fully
bidirectional bridge to python, which  transfers data buffers and does
not  copy. Many  ideas used  behind this  library have  been developed
together wit Timo Streckenbach from  the WIAS pdelib project, and they
evolved over may years. In some sense, numcxx is a rewrite of parts of
the code developed in that project.


Similar libraries
marray
eigen
armadillo
numcpp


   pybind11? Young project.
   https://community.lsst.org/t/using-pybind11-instead-of-swig-to-wrap-c-code/1096
   http://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays/  
