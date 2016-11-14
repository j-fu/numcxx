# This is a makefile in the dialect of GNU Make
# which more convenient but less general then standard
# UNIX make

WIR=/net/wir


#Which compiler to uses
CXX=g++
CXX=clang++
CXXFLAGS=-g -std=c++11 -I$(WIR)/include
LIBS=-llapack -lblas -L$(WIR)/lib -lumfpack -lamd -lcolamd -lcholmod


# implicit rules how to compile cxx files
.SUFFIXES: .cxx

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< $(LIBS) 

# List of sources and programs to be compiled
SOURCES=$(wildcard *.cxx)
PROGRAMS=$(patsubst %.cxx,%,$(SOURCES))


# Default target: compile all
all: $(PROGRAMS)

# clean-up
clean: 
	-rm $(PROGRAMS)


test:
	for f in $(PROGRAMS); do echo $$f:; ./$$f; done
