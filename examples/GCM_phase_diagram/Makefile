SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: fluid_compression_gcm
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: fluid_compression_gcm

fluid_compression_gcm: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o fluid_compression_gcm -lstdc++fs
	mv fluid_compression_gcm ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean:
	$(RM) fluid_compression_gcm *.o *.dat *.txt

clean_keep_data:
	$(RM) fluid_compression_gcm *.o