SHELL = /bin/bash
include ../../Makefile.variables
LIB = ../../lib
BIN = ../../bin
SRC = ../../src
INC = -I ../../lib -I ../../src -I ../../tools/FileIO
LN  = -L ../../lib -L ../../src -L ../../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: database_epx_pp
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: database_epx_pp

database_epx_pp: *.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) *.o $(LINKER) -o database_epx_pp -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) -o $@ $<

clean:
	$(RM) database_epx_pp *.o *.dat *.txt

clean_keep_data:
	$(RM) database_epx_pp *.o