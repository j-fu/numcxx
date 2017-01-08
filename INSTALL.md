Installation instructions
=========================

The library has been tested under
Linux
MacOSX
Win10/Cygwin

## Linux

The following packages should be installed
using the Linux package manager:

````
gcc
g++
make (GNU Make, which is the default on Linux and Mac anyway)
blas
lapack
suitesparse (umfpack)
python (2 or 3)
python-devel
python-numpy
python-matplotlib
swig
````

## Mac

We need the same packages as on Linux.

- Install Xcode from the  App-Store 

- Trigger installaion of command line developer tools in the terminal via 

````
$ gcc
````

A  dialogue   window  should  pop   up,  click  on  install

- Check with 

````
$ xcode-select -p /Library/Developer/CommandLineTools
````~

- Install [Homebrew](http://brew.sh/index.html) + [Cakebrew](https://www.cakebrew.com/)

- Install via homebrew (from science tree)
````
make
cmake
suite-sparse
````


# Cygwin

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
      make
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





