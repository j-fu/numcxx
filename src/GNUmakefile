NUMCXX=..
include $(NUMCXX)/Makefile.inc

BINARIES= libnumcxx.a _numcxxwrap.so

# default rule, makes binaries and copies them to ../lib
all: $(BINARIES)
	-mkdir -p $(NUMCXX)/lib
	-mkdir -p $(NUMCXX)/lib/python
	-mkdir -p $(NUMCXX)/lib/python/numcxx
	-cp -p libnumcxx.a $(NUMCXX)/lib
	-cp -p *.py *.so $(NUMCXX)/lib/python/numcxx

# numcxx library
libnumcxx.a:  libnumcxx.a(triangle.o) libnumcxx.a(simplegrid.o)
libnumcxx.a(triangle.o): triangle.o
libnumcxx.a(simplegrid.o): simplegrid.o


# compile triangle (C code)
triangle.o: $(NUMCXX)/triangle/triangle.c $(NUMCXX)/triangle/triangle.h
	$(CC) -DTRILIBRARY -DREAL=double -DANSI_DECLARATORS $(CFLAGS) -c  $< -o $@

# create numcxxwrap.cxx using swig
numcxxwrap.cxx: $(NUMCXX_SRC)/numcxx.i $(NUMCXX_HEADERS)
	$(SWIG)  -c++ -python  -o numcxxwrap.cxx $(NUMCXX_SRC)/numcxx.i 

# dependency for numcxxwrap.o 
numcxxwrap.o: numcxxwrap.cxx

# create shared library
_numcxxwrap.so: numcxxwrap.cxx libnumcxx.a
	$(CXX) $(CXXFLAGS) $(PYTHON_LIBS) $(LINALG_LIBS) $(SHLDFLAGS) $< libnumcxx.a -o $@

clean:
	-rm -r *.o  *.a numcxxwrap*.* *.dSYM *.pyc numcxx.py $(BINARIES) *~

