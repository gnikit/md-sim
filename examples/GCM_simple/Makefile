SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: simple_run
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: simple_run

simple_run: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o simple_run -lstdc++fs
	cp simple_run ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt

clean_keep_data:
	$(RM) simple_run *.o