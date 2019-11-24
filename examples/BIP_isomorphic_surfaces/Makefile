SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: isomorph_line
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: isomorph_line

isomorph_line: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o isomorph_line -lstdc++fs
	cp isomorph_line ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt

clean_keep_data:
	$(RM) isomorph_line *.o