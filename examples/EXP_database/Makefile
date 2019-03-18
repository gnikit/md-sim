SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: database_epx_pp
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: database_epx_pp

database_epx_pp: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o database_epx_pp -lstdc++fs
	mv database_epx_pp ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean:
	$(RM) database_epx_pp *.o *.dat *.txt

clean_keep_data:
	$(RM) database_epx_pp *.o