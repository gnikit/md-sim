SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: gcm_compression_density_increment
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: gcm_compression_density_increment

gcm_compression_density_increment: gcm_compression_density_increment.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) gcm_compression_density_increment.o $(LINKER) -o gcm_compression_density_increment -lstdc++fs
	cp gcm_compression_density_increment ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean:
	$(RM) *.o *.dat *.txt

clean_keep_data:
	$(RM) *.o