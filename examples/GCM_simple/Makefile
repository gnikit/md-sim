SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: gcm_simple_run
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: gcm_simple_run

gcm_simple_run: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o gcm_simple_run -lstdc++fs
	cp gcm_simple_run ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt *.log

clean_keep_data:
	$(RM) gcm_simple_run *.o