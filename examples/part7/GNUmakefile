


# This needs to be customized to point to the numcxx installation
# on your system
NUMCXX=../..


ifeq ($(wildcard $(NUMCXX)/Makefile.inc),)
NUMCXX=/net/wir/numcxx
endif


include $(NUMCXX)/Makefile.inc
CXXFLAGS:=$(CXXFLAGS) -fopenmp -pthread
CFLAGS:=$(CFLAGS) -pthread



PROGRAMS=\
	pthreads\
	cpp11threads\
	cpp11mutex\
	cppomp\
	cppomp-critical\
	omp-scalpro-wrong\
	omp-scalpro-critical\
	omp-scalpro-reduction\
	mpi-sum\
	mpi-hello

#
# Special targets for mpi parallel stuff
# 

mpi-hello:	mpi-hello.cxx
	mpic++ -std=c++11 -O3 $< -o $@

mpi-sum:	mpi-sum.cxx
	mpic++ -std=c++11 -O3 $< -o $@


BINARIES= $(PROGRAMS) 
all: $(BINARIES)




clean:
	-rm -r *.o  swigwrap*.* *.dSYM *.pyc __pycache__ convdiff.py $(BINARIES) *~
