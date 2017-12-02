NUMCXX Installation
===================

The library has been tested under
Linux,
MacOSX,
Win10/Cygwin

## Prerequisites

Installation of the python components is currently not necessary.

### Linux

The following packages should be installed
using the Linux package manager (from system to system,
they have different names)

````
gcc
g++
cmake 
blas
lapack 
suitesparse 
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

  - Install the cygwin packages

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

## Installation of numcxx

### TU Berlin UNIX Pool

Here,  the  current  version  of  numcxx,  without  the  examples,  is
installed in the directory /net/wir  The example subdirectories can be
copied separately as a whole.

### Virtual Machine

The debian-numcxx  virtual machine can  be downloaded from  the course
home page. Load this machine into VirtualBox and start it.
Vtk, numcxx and vtkfig and codeblocks are installed on this machine. Examples
can be run with the numcxx-build utility.

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

### Working with the examples


First, copy an example directory to your working directory

````
$ cp -r <numcxx>/examples/whatever .
$ cd whatever
````

Each of the example subdirectories has its own CMakeLists.txt which
looks for an installed numcxx during the setup phase.

### Plain unix


Set up a build subdirectory and compile the project:

````
$ mkdir .build
$ cd .build
$ CXXFLAGS=-std=c++11 cmake ..
$ cmake --build .
$ ./<whatever example>
````

### Code::Blocks

Set up a build subdirectory for codeblocks:
````
$ mkdir .build
$ cd .build
$ CXXFLAGS=-std=c++11 cmake -G"CodeBlocks - Unix Makefiles" ..
$ codeblocks <whatver the project name is>.cbp
````


In codeblocks, select a target under "Build/Select Target".
In the management pane (left), select the corresponding
source file. Now you can build and run   via the correspondin
menu items of codeblocks.

### Other IDEs
Many of them should work. Like with  code blocks, find out the name of
the generator via ``cmake -G``. If you don't find your IDE name there,
there might be  some chance that your IDE supports  cmake directly (as
e.g. QTCreator does).

