# This needs to be customized to point to the numcxx installation
# on your system
NUMCXX=../..


ifeq ($(wildcard $(NUMCXX)/Makefile.inc),)
NUMCXX=/net/wir/numcxx
endif


include $(NUMCXX)/Makefile.inc

PROGRAMS=
BINARIES= $(PROGRAMS) _convdiff.$(SO) _fvm2d.$(SO)
all: $(BINARIES)


swigwrap_convdiff.cxx: convdiff.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_convdiff.cxx convdiff.i

_convdiff.$(SO): swigwrap_convdiff.cxx convdiff.hxx
	$(CXX) $(CXXFLAGS) $< $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) -o $@


swigwrap_fvm2d.cxx: fvm2d.i fvm2d.hxx $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o swigwrap_fvm2d.cxx fvm2d.i

_fvm2d.$(SO): swigwrap_fvm2d.cxx fvm2d.hxx
	$(CXX) $(CXXFLAGS) $< $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) -o $@




clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc __pycache__ convdiff.py $(BINARIES) *~
