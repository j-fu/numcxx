CXX=g++-5
#CXX=clang++
PYTHON=python
CXXFLAGS=-fPIC -O3 -g -DNUMCXX_CHECK_BOUNDS -std=c++11 `$(PYTHON) find_python_paths.py --includes`
PYTHONLIBS=`$(PYTHON) find_python_paths.py --libs`
LALIBS=-llapack -lblas
SHLDFLAGS= -shared 

.SUFFIXES: .cxx

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LALIBS)

all: test-numcxx _numcxx.so


numcxx.o:  numcxx.hxx
numcxx.o:  numcxx.cxx


numcxx_wrap.o: numcxx_wrap.cxx

test-numcxx: test-numcxx.cxx numcxx.hxx

numcxx_wrap.cxx: numcxx.i numcxx.hxx
	swig  -c++ -python   numcxx.i 

_numcxx.so: numcxx_wrap.o 
	$(CXX) $(CXXFLAGS) $(PYTHONLIBS) $(SHLDFLAGS) $< -o $@

