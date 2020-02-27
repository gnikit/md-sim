SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: gcm_solid_data_generator
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: gcm_solid_data_generator

gcm_solid_data_generator: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o gcm_solid_data_generator -lstdc++fs
	cp gcm_solid_data_generator ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt *.log

clean_keep_data:
	$(RM) gcm_solid_data_generator *.o