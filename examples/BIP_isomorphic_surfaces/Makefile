SHELL = /bin/bash
include ../../Makefile.variables
LIB = ../../lib
BIN = ../../bin
SRC = ../../src
INC = -I ../../lib -I ../../src -I ../../tools/FileIO
LN  = -L ../../lib -L ../../src -L ../../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: isomorph_line
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: isomorph_line

isomorph_line: *.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) *.o $(LINKER) -o isomorph_line -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) -o $@ $<

clean:
	$(RM) isomorph_line *.o *.dat *.txt

clean_keep_data:
	$(RM) isomorph_line *.o