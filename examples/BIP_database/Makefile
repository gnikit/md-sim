SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: bip_database
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: bip_database

bip_database: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o bip_database -lstdc++fs

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt *.log

clean_keep_data:
	$(RM) bip_database *.o