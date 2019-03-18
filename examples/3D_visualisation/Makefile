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
debug: visualise_fluid

visualise_fluid: *.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) *.o $(LINKER) -o visualise_fluid -lstdc++fs
	mv visualise_fluid ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean:
	$(RM) visualise_fluid *.o *.dat *.txt

clean_keep_data:
	$(RM) visualise_fluid *.o