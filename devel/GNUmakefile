NUMCXX=..
include $(NUMCXX)/Makefile.inc

PROGRAMS=test-sparse test-et  test-numcxx3  featuretest test-numcxx

BINARIES= $(PROGRAMS) 
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

test-triangle: test-triangle.cxx $(NUMCXX_HEADERS) triangle.o
	$(CXX) $(CXXFLAGS) -o $@  $< triangle.o $(LINALG_LIBS)

test-numcxx: test-numcxx.cxx $(NUMCXX_HEADERS) test-numcxx2.o
	$(CXX) $(CXXFLAGS) -o $@  $< test-numcxx2.o $(LINALG_LIBS)
clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc numcxx.py $(BINARIES) *~


