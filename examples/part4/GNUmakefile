NUMCXX=../..
include $(NUMCXX)/Makefile.inc

BINARIES= $(PROGRAMS) _mymodule.so
all: $(BINARIES)


swigwrap_mymodule.cxx: mymodule.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_mymodule.cxx mymodule.i

_mymodule.so: swigwrap_mymodule.cxx 
	$(CXX) $(CXXFLAGS) $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) $< -o $@
test: 
	PYTHONPATH=$(NUMCXX)/lib/python $(PYTHON) test.py

clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc __pycache__ mymodule.py $(BINARIES) *~
