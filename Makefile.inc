###########################################
# Make variables to adapt. Consider to put
# any changes to MyMakefile.inc

# C++ compiler
CXX=g++

# C compiler
CC=gcc

# Python executable
PYTHON=python

# C++ compiler flags
CXXFLAGS=-O3  -std=c++11 

# Control array bounds check in numcxx
CHECKBOUNDS_FLAGS=-DNUMCXX_CHECK_BOUNDS

# C compiler flags
CFLAGS=-O3 

# Compiler flags for generating position independent code
# (for Mac/Unix shared objects, ignored on windows)
PICFLAGS=-fPIC

# Linker flag for creating shared libraries
SHLDFLAGS= -shared 

# Dense+Sparse linear algebra libraries
LINALG_LIBS=-L/usr/local/lib -llapack -lblas -lumfpack -lamd -lcolamd -lcholmod -lsuitesparseconfig

# Sparse linear algebra includes
LINALG_INCLUDES=-I/usr/include/suitesparse -I/usr/local/include 


# Possible locations for swig
_SWIG1=/usr/bin/swig
_SWIG2=/usr/bin/swig3.0
_SWIG3=/usr/local/bin/swig
_SWIG4=/Users/fuhrmann/Wias/work/pdelib2/external/lib.gnu_hb/swig-3.0.5-2/bin/swig
_SWIG5=/Users/fuhrmann/Wias/work/pdelib2/external/lib.gnu_hb/swig-3.0.5-2/bin/swig


# End of customizable part
###############################################################

###############################################################
# This file is meant to be included from subdirectory makefiles
# which define $(NUMCXX). It uses GNUMake syntax (which 
# defines ':=' and 'ifeq' in addition  to the usual makefile 
# syntax.

# If exists, include a file  where you can overwrite 
# any of the previous flags without changing
# the present file file.
sinclude $(NUMCXX)/MyMakefile.inc

# Include directory for numcxx 
NUMCXX_INCLUDE=$(NUMCXX)/include

# Source directory for numcxx
NUMCXX_SRC=$(NUMCXX)/src

CXXFLAGS:=$(CXXFLAGS) $(CHECKBOUNDS_FLAGS)

# numcxx header files
NUMCXX_HEADERS=$(wildcard $(NUMCXX_INCLUDE)/*/*.hxx $(NUMCXX_INCLUDE)/*/*.ixx)

# Amend cxxflags by necessary includes
CXXFLAGS:=$(CXXFLAGS)   -I$(NUMCXX)/include $(LINALG_INCLUDES) `$(PYTHON) $(NUMCXX_SRC)/find_python_paths.py --includes`

# Amend libs by necessary python libs
PYTHON_LIBS=`$(PYTHON) $(NUMCXX_SRC)/find_python_paths.py --libs`


# Default rules
.SUFFIXES: .c .cxx .cpp

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(NUMCXX)/lib/libnumcxx.a $(LINALG_LIBS) 


SWIG="error: swig_does_not_exist!"

# Figure out swig executable
ifneq ($(wildcard $(_SWIG1)),)
SWIG=$(_SWIG1)
endif

ifneq ($(wildcard $(_SWIG2)),)
SWIG=$(_SWIG2)
endif

ifneq ($(wildcard $(_SWIG3)),)
SWIG=$(_SWIG3)
endif

ifneq ($(wildcard $(_SWIG4)),)
SWIG=$(_SWIG4)
endif

ifneq ($(wildcard $(_SWIG5)),)
SWIG=$(_SWIG5)
endif



# Figure out extension of shared objects
ifneq ($(wildcard /usr/bin/cygwin*.dll),)
SO=dll
else
SO=so
CXXFLAGS:=$(CXXFLAGS) $(PICFLAGS)
CFLAGS:=$(CFLAGS) $(PICFLAGS)
endif
