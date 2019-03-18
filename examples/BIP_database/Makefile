SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: database
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: database

database: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o database -lstdc++fs
	mv database ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean:
	$(RM) database *.o *.dat *.txt

clean_keep_data:
	$(RM) database *.o