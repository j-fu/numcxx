# This needs to be customized to point to the numcxx installation
# on yoru system
NUMCXX=../..

ifeq ($(wildcard $(NUMCXX)/Makefile.inc),)
NUMCXX=/net/wir	
endif

include $(NUMCXX)/Makefile.inc

BINARIES= $(PROGRAMS) _mymodule.$(SO)
all: $(BINARIES)


swigwrap_mymodule.cxx: mymodule.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_mymodule.cxx mymodule.i

_mymodule.$(SO): swigwrap_mymodule.cxx 
	$(CXX) $(CXXFLAGS) $< $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) -o $@
test: 
	PYTHONPATH=$(NUMCXX)/lib/python $(PYTHON) featuretest.py

clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc __pycache__ mymodule.py $(BINARIES) *~
