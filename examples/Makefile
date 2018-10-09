SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
SRC = ../src
INC = -I ../lib -I ../src -I ../tools/FileIO
LN = -L ../lib -L ../src -L ../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: database isomorph fluid_compression $(LIB)/libmd.a

database: database.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) database.o $(LINKER) -o database
	
isomorph: isomorph_line.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) isomorph_line.o $(LINKER) -o isomorph

fluid_compression: fluid_compression.o $(LIB)/libmd.a 
	$(CXX) $(FLAGS) fluid_compression.o $(LINKER) -o fluid_compression

database.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) make_database.cpp -o database.o

isomorph_line.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) isomorph_line.cpp -o isomorph_line.o

fluid_compression.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) fluid_compression.cpp -o fluid_compression.o


clean:
	$(RM) isomorph database fluid_compression *.o 
	$(RM) ./example_data/*.txt
	$(RM) ../data/initial_velocities*
	$(RM) ../data/vx* ../data/vy* ../data/vz*


clean_keep_data:
	$(RM) isomorph database fluid_compression *.o
