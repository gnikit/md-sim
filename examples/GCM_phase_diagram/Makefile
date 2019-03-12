SHELL = /bin/bash
include ../../Makefile.variables
LIB = ../../lib
BIN = ../../bin
SRC = ../../src
INC = -I ../../lib -I ../../src -I ../../tools/FileIO
LN  = -L ../../lib -L ../../src -L ../../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: fluid_compression
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: fluid_compression

fluid_compression: *.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) *.o $(LINKER) -o fluid_compression -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) -o $@ $<

clean:
	$(RM) fluid_compression *.o *.dat *.txt

clean_keep_data:
	$(RM) fluid_compression *.o