# This is a makefile in the dialect of GNU Make
# which more convenient but less general then standard
# UNIX make

#Which compiler to uses
CXX=g++
CXXFLAGS=-g -std=c++11

# implcit rules how to compile cxx files
.SUFFIXES: .cxx

.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< 

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

