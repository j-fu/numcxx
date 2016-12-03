CXX=g++-5
CC=gcc-5
CXX=clang++
CC=clang
CXX=g++
CC=gcc
#CXX=icpc
PYTHON=python

NUMCXX=../
NUMCXX_INCLUDE=$(NUMCXX)/include
NUMCXX_PYTOOLS=$(NUMCXX)/pytools

CFLAGS=-fPIC -O0 -g -DTRILIBRARY -DREAL=double -DANSI_DECLARATORS 

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

PROGRAMS=test-sparse test-et  test-numcxx3  featuretest test-numcxx

BINARIES= $(PROGRAMS) _numcxx.so 
all: $(BINARIES)

test:
	for f in $(PROGRAMS); do echo $$f:; ./$$f; done


NUMCXX_HEADERS=$(wildcard $(NUMCXX_INCLUDE)/*/*.hxx $(NUMCXX_INCLUDE)/*/*.ixx)


test-numcxx3: $(NUMCXX_HEADERS)
test-numcxx2.o: $(NUMCXX_HEADERS)
featuretest: $(NUMCXX_HEADERS)
test-et: $(NUMCXX_HEADERS)
test-sparse: $(NUMCXX_HEADERS)
triangle.o: $(NUMCXX)/triangle/triangle.c

swigwrap_numcxx.o: swigwrap_numcxx.cxx

test-triangle: test-triangle.cxx $(NUMCXX_HEADERS) triangle.o
	$(CXX) $(CXXFLAGS) -o $@  $< triangle.o $(LINALG_LIBS)

test-numcxx: test-numcxx.cxx $(NUMCXX_HEADERS) test-numcxx2.o
	$(CXX) $(CXXFLAGS) -o $@  $< test-numcxx2.o $(LINALG_LIBS)

swigwrap_numcxx.cxx: $(NUMCXX_PYTOOLS)/numcxx.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_numcxx.cxx $(NUMCXX_PYTOOLS)/numcxx.i 

_numcxx.so: swigwrap_numcxx.cxx triangle.o
	$(CXX) $(CXXFLAGS) $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) $< triangle.o -o $@

clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc numcxx.py $(BINARIES) *~


