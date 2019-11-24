SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: gcm_compression
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: gcm_compression

gcm_compression: gcm_compression.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) gcm_compression.o $(LINKER) -o gcm_compression -lstdc++fs
	cp gcm_compression ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.o *.dat *.txt 

clean_keep_data:
	$(RM) *.o gcm_compression