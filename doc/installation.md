NUMCXX Installation
===================

The library has been tested under
Linux,
MacOSX,
Win10/Cygwin

## Prerequisites

Installation of the python components is only necessary when
the experimental python interface is needed.

### Linux

The following packages should be installed
using the Linux package manager (from system to system,
they have different names)

Mandatory:
````
gcc
g++
cmake 
blas
lapack 
suitesparse 
````

Optional:
````
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

- A  dialogue   window  should  pop   up,  click  on  `install`

````
$ xcode-select -p /Library/Developer/CommandLineTools
````

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

Mandatory:

````
      gcc-g++
      wget
      curl
      openblas
      cmake
      lapack
      libumfpack-devel
````


- If you get  tired clicking into the setup gui, download [apt-cyg](https://raw.githubusercontent.com/transcode-open/apt-cyg/master/apt-cyg)
  and copy the script to ``/usr/bin`` (which is ``C:\cygwin64\usr\bin`` from the
  DOS prompt). Then you can install all the packages using ``apt-cyg install``.

- This should be sufficient to run make and make test.
  Before that, close the cygwin shell and open it again in order
  to catch the proper path settings for ``lapack``

- For Python graphics (matplotlib), proceed further:

  - Install the cygwin packages

````
  python
  python-numpy
  python-devel 
  swig
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

## Installation of numcxx

### Build on your own system

Go to the  numcxx root directory and set up the cmake project:

````
$ mkdir .build
$ cd .build
$ CXXFLAGS=-std=c++11 cmake ..
$ cmake --build .
````

This should compile the library and make it findeable by depending
projects via an entry im the user CMake project registry found
under ``.cmake`` in your home directory.

Moreover, it places the ``numcxx-build`` script into the subdirectory ``bin``.
This script can be used to build the examples and small projects.


### TU Berlin UNIX Pool

Here,  the  current  version  of  numcxx  is
installed in the directory ``/net/wir``

### Virtual Machine

The debian-numcxx  virtual machine can  be downloaded from  the course
home page. Load this machine into VirtualBox and start it.
Vtk, numcxx and vtkfig and codeblocks are installed on this machine. Examples
can be run with the numcxx-build utility.


### Working with the examples


First, copy an example directory to your working directory

````
$ cp -r <numcxx>/examples/whatever .
$ cd whatever
````

All the examples can be built with the numcxx-build script.


