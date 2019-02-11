SHELL = /bin/bash
include ../../Makefile.variables
LIB = ../../lib
BIN = ../../bin
SRC = ../../src
INC = -I ../../lib -I ../../src -I ../../tools/FileIO
LN  = -L ../../lib -L ../../src -L ../../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: visualise_fluid
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: visualise_fluid

visualise_fluid: *.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) *.o $(LINKER) -o visualise_fluid -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) -o $@ $<

clean:
	$(RM) visualise_fluid *.o *.dat *.txt

clean_keep_data:
	$(RM) visualise_fluid *.o