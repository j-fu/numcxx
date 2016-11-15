CXX=g++-5
#CXX=g++
CXX=clang++
#CXX=icpc
PYTHON=python

NUMCXX=../
NUMCXX_INCLUDE=$(NUMCXX)/include
NUMCXX_PYTOOLS=$(NUMCXX)/pytools

CXXFLAGS=-fPIC -O3 -g  -DNUMCXX_CHECK_BOUNDS -I$(NUMCXX)/include -I/usr/include/suitesparse -I/usr/local/include -std=c++11 `$(PYTHON) $(NUMCXX_PYTOOLS)/find_python_paths.py --includes`
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

swigwrap_numcxx.o: swigwrap_numcxx.cxx

test-numcxx: test-numcxx.cxx $(NUMCXX_HEADERS) test-numcxx2.o
	$(CXX) $(CXXFLAGS) -o $@  $< test-numcxx2.o $(LINALG_LIBS)

swigwrap_numcxx.cxx: $(NUMCXX_PYTOOLS)/numcxx.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_numcxx.cxx $(NUMCXX_PYTOOLS)/numcxx.i 

_numcxx.so: swigwrap_numcxx.cxx
	$(CXX) $(CXXFLAGS) $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) $< -o $@

clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc numcxx.py $(BINARIES) *~


