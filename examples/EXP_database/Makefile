SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: exp_database
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: exp_database

exp_database: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o exp_database -lstdc++fs

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt

clean_keep_data:
	$(RM) exp_database *.o