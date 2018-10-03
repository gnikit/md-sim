SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
SRC = ../src
# TODO: This will obviously not work if isomorph.o, MD.o stat_analysis.o do not exist
# 			at this point 'make all' has to be called and 'database' and 'iso' cannot be called independently
all: make_database isomorphic_surface fluid_compression $(LIB)/libmd.a

make_database: make_database.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) make_database.o $(LINKER) -o make_database
	
isomorphic_surface: isomorph_line.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) isomorph_line.o $(LINKER) -o isomorphic_surface

fluid_compression: fluid_compression.o $(LIB)/libmd.a 
	$(CXX) $(FLAGS) fluid_compression.o $(LINKER) -o fluid_compression

make_database.o: make_database.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) make_database.cpp -o make_database.o

isomorph_line.o: isomorph_line.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) isomorph_line.cpp -o isomorph_line.o

fluid_compression.o: fluid_compression.cpp
	$(CXX) $(CXXFLAGS) -I$(SRC) fluid_compression.cpp -o fluid_compression.o


clean:
	$(RM) isomorphic_surface make_database fluid_compression *.o 
	$(RM) ./example_data/*.txt

clean_keep_data:
	$(RM) isomorphic_surface make_database fluid_compression *.o
