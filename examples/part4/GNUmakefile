CXX=g++
CC=gcc
PYTHON=python3

NUMCXX=../../
NUMCXX_INCLUDE=$(NUMCXX)/include
NUMCXX_PYTOOLS=$(NUMCXX)/pytools

CFLAGS=-fPIC -O3 -g -DTRILIBRARY -DREAL=double -DANSI_DECLARATORS 
CXXFLAGS=$(CFLAGS)  -DNUMCXX_CHECK_BOUNDS -I$(NUMCXX)/include -I/usr/include/suitesparse -I/usr/local/include -std=c++11 `$(PYTHON) $(NUMCXX_PYTOOLS)/find_python_paths.py --includes`
PYTHON_LIBS=`$(PYTHON) $(NUMCXX_PYTOOLS)/find_python_paths.py --libs`
LINALG_LIBS=-L/usr/local/lib -llapack -lblas -lumfpack -lamd -lcolamd -lcholmod -lsuitesparseconfig


# try for swig
_SWIG1=/usr/bin/swig
_SWIG2=/usr/bin/swig3.0
_SWIG3=/Users/fuhrmann/Wias/work/pdelib2/external/lib.gnu_hb/swig-3.0.5-2/bin/swig


ifneq ($(wildcard $(_SWIG1)),)
SWIG=$(_SWIG1)
endif

ifneq ($(wildcard $(_SWIG2)),)
SWIG=$(_SWIG2)
endif

ifneq ($(wildcard $(_SWIG3)),)
SWIG=$(_SWIG3)
endif


SHLDFLAGS= -shared 

.SUFFIXES: .cxx


.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LINALG_LIBS)

PROGRAMS=

BINARIES= $(PROGRAMS) _mymodule.so
all: $(BINARIES)


NUMCXX_HEADERS=$(wildcard $(NUMCXX_INCLUDE)/*/*.hxx $(NUMCXX_INCLUDE)/*/*.ixx)



swigwrap_mymodule.cxx: mymodule.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_mymodule.cxx mymodule.i

_mymodule.so: swigwrap_mymodule.cxx 
	$(CXX) $(CXXFLAGS) $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) $< -o $@


clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc __pycache__ mymodule.py $(BINARIES) *~



