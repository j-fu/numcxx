CXX=g++-5
CXX=g++
CXX=clang++
#CXX=icpc
NUMCXX=../include/numcxx

PYTHON=python
CXXFLAGS=-fPIC -O3 -g  -DNUMCXX_CHECK_BOUNDS -I ../include -std=c++11 `$(PYTHON) $(NUMCXX)/find_python_paths.py --includes`
PYTHONLIBS=`$(PYTHON) $(NUMCXX)/find_python_paths.py --libs`
LALIBS=-llapack -lblas
SHLDFLAGS= -shared 

.SUFFIXES: .cxx

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LALIBS)


BINARIES=test-sparse test-et  test-numcxx3  featuretest test-numcxx _numcxx.so 
all: $(BINARIES)


NUMCXX_HEADERS=$(wildcard $(NUMCXX)/*.hxx)


test-numcxx3: $(NUMCXX_HEADERS)
test-numcxx2: $(NUMCXX_HEADERS)
featuretest: $(NUMCXX_HEADERS)
test-et: $(NUMCXX_HEADERS)
test-sparse: $(NUMCXX_HEADERS)

swigwrap_numcxx.o: swigwrap_numcxx.cxx

test-numcxx: test-numcxx.cxx $(NUMCXX_HEADERS) test-numcxx2.o
	$(CXX) $(CXXFLAGS) -o $@  $< test-numcxx2.o $(LALIBS)

swigwrap_numcxx.cxx: $(NUMCXX)/numcxx.i $(NUMCXX_HEADERS)
	swig  -c++ -python  -o swigwrap_numcxx.cxx $(NUMCXX)/numcxx.i 

_numcxx.so: swigwrap_numcxx.cxx
	$(CXX) $(CXXFLAGS) $(PYTHONLIBS) $(SHLDFLAGS) $< -o $@

clean:
	-rm *.o  swigwrap*.*  *.pyc numcxx.py $(BINARIES) *~

