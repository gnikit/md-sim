SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: exp_simple_run
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: exp_simple_run

exp_simple_run: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o exp_simple_run -lstdc++fs
	cp exp_simple_run ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt

clean_keep_data:
	$(RM) exp_simple_run *.o