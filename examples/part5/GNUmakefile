NUMCXX=../..
include $(NUMCXX)/Makefile.inc

BINARIES= $(PROGRAMS) _fem2d.$(SO)
all: $(BINARIES)


swigwrap_fem2d.cxx: fem2d.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_fem2d.cxx fem2d.i

_fem2d.$(SO): swigwrap_fem2d.cxx fem2d.hxx
	$(CXX) $(CXXFLAGS) $< $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) -o $@



clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc __pycache__ fem2d.py $(BINARIES) *~
