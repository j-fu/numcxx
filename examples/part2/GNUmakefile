# This is a Makefile in the dialect of GNU Make
# which is more convenient but less general then standard
# UNIX make

## Make variables
# Location of installation in TUB unix pool
WIR=/net/wir

# Location within numcxx distribution
NUMCXX=../..


#Which compiler to use
CXX=g++

#Other choices: 
#clang++
#g++-5

# Include directories for numcxx library
NUMCXX_INCLUDES=-I$(NUMCXX)/include -I$(WIR)/numcxx/include
#                ^ numcxx distribution
#                                ^ UNIX pool

# Include directories for suitesparse library
SUITESPARSE_INCLUDES=-I/usr/include/suitesparse -I/usr/local/include

#Flags to be used during compilation
CXXFLAGS=-g -O3 -std=c++11 $(NUMCXX_INCLUDES) $(SUITESPARSE_INCLUDES)
#         ^ Generate code for debugging
#            ^ optimize code
#               ^ use the C++11 standard



# libraries to be linked with
# indicating -lname triggers a search for
# -libname.a or -libname.so on system search paths
LIBS=-llapack -lblas
#      ^ lapack Linear Algebar PACKage  containing tools for LU factorization of dense matrices etc.
#               ^ Basic Linear Subroutines PACKage containing vector and matrix x vector operations

## Implicit rules
# implcit rules how to compile cxx files
.SUFFIXES: .cxx

# How to compile object file from cxx file
.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@

# How to compile and link program
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LIBS) 
#                                     ^ Link with the libraries indicated in libs
## List of sources and programs
# List of sources  to be compiled
SOURCES=$(wildcard *.cxx)

# List of programs to be generated
PROGRAMS=$(patsubst %.cxx,%,$(SOURCES))

## Targets
# Default target: compile all programs
# Here, make checks if  a source is newer than a program
# for all programs in the list.
# In this case, it tries to find an implicit rule to 
# generate the source from the program and applies it.
# Here, it will be ".cxx".
all: $(PROGRAMS)



# Remove all programs
clean: 
	-rm $(PROGRAMS)

#run  all programs
test:
	for f in $(PROGRAMS); do echo $$f:; ./$$f; done

