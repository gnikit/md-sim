SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: fluid_compression_bip
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: fluid_compression_bip

fluid_compression_bip: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o fluid_compression_bip -lstdc++fs
	mv fluid_compression_bip ../examplebin/


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean:
	$(RM) fluid_compression_bip *.o *.dat *.txt

clean_keep_data:
	$(RM) fluid_compression_bip *.o