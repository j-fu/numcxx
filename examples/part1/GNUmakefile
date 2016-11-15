# This is a Makefile in the dialect of GNU Make
# which more convenient but less general then standard
# UNIX make

## Make variables
#Which compiler to use
CXX=g++

#Other choices: 
#clang++
#g++-5

#Flags to be used during compilation
CXXFLAGS=-g -O3 -std=c++11
#         ^ Generate code for debugging
#            ^ optimize code
#               ^ use the C++11 standard


## Implicit rules
# implcit rules how to compile cxx files
.SUFFIXES: .cxx

# How to compile object file from cxx file
.cxx.o:
	$(CXX) $(CXXFLAGS) -c  $< -o $@

# How to compile and link program
.cxx:
	$(CXX) $(CXXFLAGS) -o $@  $< 

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

