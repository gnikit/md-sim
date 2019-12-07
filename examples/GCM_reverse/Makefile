SHELL = /bin/bash
include ../../Makefile.variables
LIB := ../../lib
BIN := ../../bin
SRC := ../../src
INC := -I ../../include
LN  := -L ../../lib

LINKER = $(INC) $(LN) -lmd

all: gcm_reverse
			
debug: CXXFLAGS = $(DCXXFLAGS)
debug: gcm_reverse

gcm_reverse: gcm_reverse.o $(LIB)/libmd.a
	$(CXX) $(CXXFLAGS) gcm_reverse.o $(LINKER) -o gcm_reverse -lstdc++fs
	cp gcm_reverse ../examplebin/

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SRC) $<

clean: clean_keep_data
	$(RM) *.dat *.txt *.log 

clean_keep_data:
	$(RM) *.o gcm_reverse