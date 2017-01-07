###########################################
# Make variables to adapt

# C++ compiler
CXX=g++

# C compiler
CC=gcc

# Python executable
PYTHON=python3

# C++ compiler flags
CXXFLAGS=-fPIC -O3  -std=c++11 

#Check bounds flag for numcxx
CXXFLAGS:=$(CXXFLAGS) -DNUMCXX_CHECK_BOUNDS

# C compiler flags
CFLAGS=-fPIC -O3 

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


# end of customizable part
###############################################################
# This file is meant to be included from subdirectory makefiles
# which define $(NUMCXX)


# include directory for numcxx 
NUMCXX_INCLUDE=$(NUMCXX)/include

# source directory for numcxx
NUMCXX_SRC=$(NUMCXX)/src

# numcxx header files
NUMCXX_HEADERS=$(wildcard $(NUMCXX_INCLUDE)/*/*.hxx $(NUMCXX_INCLUDE)/*/*.ixx)

# amend cxxflags by necessary includes
CXXFLAGS:=$(CXXFLAGS)   -I$(NUMCXX)/include $(LINALG_INCLUDES) `$(PYTHON) $(NUMCXX_SRC)/find_python_paths.py --includes`

# amend libs by necessary python libs
PYTHON_LIBS=`$(PYTHON) $(NUMCXX_SRC)/find_python_paths.py --libs`


# default rules
.SUFFIXES: .c .cxx .cpp

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LINALG_LIBS)


# figure out swig executable
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
