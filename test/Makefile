CXX=g++-5
#CXX=g++
CXX=clang++
#CXX=icpc
PYTHON=python

NUMCXX=../
NUMCXX_INCLUDE=$(NUMCXX)/include
NUMCXX_PYTOOLS=$(NUMCXX)/pytools

CXXFLAGS=-fPIC -O3 -g  -DNUMCXX_CHECK_BOUNDS -I ../include -std=c++11 `$(PYTHON) $(NUMCXX_PYTOOLS)/find_python_paths.py --includes`
PYTHON_LIBS=`$(PYTHON) $(NUMCXX_PYTOOLS)/find_python_paths.py --libs`
LINALG_LIBS=-llapack -lblas -lumfpack -lamd -lcolamd -lcholmod


SHLDFLAGS= -shared 

.SUFFIXES: .cxx

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LINALG_LIBS)


BINARIES=test-sparse test-et  test-numcxx3  featuretest test-numcxx _numcxx.so 
all: $(BINARIES)


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
	swig  -c++ -python  -o swigwrap_numcxx.cxx $(NUMCXX_PYTOOLS)/numcxx.i 

_numcxx.so: swigwrap_numcxx.cxx
	$(CXX) $(CXXFLAGS) $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) $< -o $@

clean:
	-rm *.o  swigwrap*.*  *.pyc numcxx.py $(BINARIES) *~


