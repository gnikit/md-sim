SHELL = /bin/bash
include ../../Makefile.variables
LIB = ../../lib
BIN = ../../bin
SRC = ../../src
INC = -I ../../lib -I ../../src -I ../../tools/FileIO
LN  = -L ../../lib -L ../../src -L ../../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: simple_run
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: simple_run

simple_run: *.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) *.o $(LINKER) -o simple_run -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) -o $@ $<

clean:
	$(RM) simple_run *.o *.dat *.txt

clean_keep_data:
	$(RM) simple_run *.o