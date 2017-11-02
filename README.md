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
  - triangle

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

## Installation 

The library has been tested under
Linux,
MacOSX,
Win10/Cygwin

Installation of the python components is currently not necessary.

### Linux

The following packages should be installed
using the Linux package manager:

````
gcc
g++
cmake 
blas
lapack
suitesparse (umfpack)
python
python-devel
python-numpy
python-matplotlib
swig
````

### Mac

We need the same packages as on Linux.

- Install Xcode from the  App-Store 

- Trigger installaion of command line developer tools in the terminal via 

````
$ gcc
````

- Check with 

A  dialogue   window  should  pop   up,  click  on  `install`

````
$ xcode-select -p /Library/Developer/CommandLineTools
````~

- Install [Homebrew](http://brew.sh/index.html) + [Cakebrew](https://www.cakebrew.com/)

- Install via homebrew (from science tree)
````
cmake
suite-sparse
````


### Cygwin

Thanks to the autor of [this post](https://bynario.com/2016-10-01-minimal-cygwin-python-data-science-installation.html)
for showing the way for matplotlib.


- Install cygwin (basic system) following
  the usual setup instructions given on [cygwin](https://www.cygwin.com/)

- Using setup_x86_64.exe, install the following additonal packages:

````
      gcc-g++
      wget
      curl
      openblas
      cmake
      lapack
      libumfpack-devel

      python
      python-numpy
      python-devel 
      swig
````



- If you get  tired clicking into the setup gui, download [apt-cyg](https://raw.githubusercontent.com/transcode-open/apt-cyg/master/apt-cyg)
  and copy the script to ``/usr/bin`` (which is ``C:\cygwin64\usr\bin`` from the
  DOS prompt). Then you can install all the packages using ``apt-cyg install``.

- This should be sufficient to run make and make test.
  Before that, close the cygwin shell and open it again in order
  to catch the proper path settings for ``lapack``

- For Python graphics (matplotlib), proceed further:

     - install the cygwin packages

    ````
      libfreetype-devel 
      python-pyqt4
      xorg-server
      xinit
    ````

     - Obtain pip, the python package manager

    ````
      $ wget https://bootstrap.pypa.io/get-pip.py
      $ python get-pip.py
    ```` 

     - Install matplotlib via pip:

    ````
      $ pip install matplotlib
    ````

    -  Start the X server and set the DISPLAY environment variable

    ````
     $ startxwin
     $ export DISPLAY=:0.0
    ````

    - use qt4agg as python backend:

    ````
     $ export MPLBACKEND=qt4agg
    ````

