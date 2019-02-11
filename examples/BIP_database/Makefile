SHELL = /bin/bash
include ../../Makefile.variables
LIB = ../../lib
BIN = ../../bin
SRC = ../../src
INC = -I ../../lib -I ../../src -I ../../tools/FileIO
LN  = -L ../../lib -L ../../src -L ../../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: database
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: database

database: *.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) *.o $(LINKER) -o database -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) -o $@ $<

clean:
	$(RM) database *.o *.dat *.txt

clean_keep_data:
	$(RM) database *.o