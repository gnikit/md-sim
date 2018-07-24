SHELL = /bin/bash
include ../Makefile.variables
LIB =../lib
BIN =../bin
# TODO: This will obviously not work if isomorph.o, MD.o stat_analysis.o do not exist
# 			at this point 'make all' has to be called and 'database' and 'iso' cannot be called independently
all: database isomorphic_surface fluid_compression $(LIB)/libmd.a

database: make_database.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) make_database.o $(LINKER) -o database

isomorphic_surface: isomorph_line.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) isomorph_line.o $(LINKER) -o isomorphic_surface

fluid_compression: fluid_compression.o $(LIB)/libmd.a 
	$(CXX) $(FLAGS) fluid_compression.o $(LINKER) -o fluid_compression

make_database.o: make_database.cpp
	$(CXX) $(CXXFLAGS) -I$(LIB) make_database.cpp -o make_database.o

isomorph_line.o: isomorph_line.cpp
	$(CXX) $(CXXFLAGS) -I$(LIB) isomorph_line.cpp -o isomorph_line.o

fluid_compression.o: fluid_compression.cpp
	$(CXX) $(CXXFLAGS) -I$(LIB) fluid_compression.cpp -o fluid_compression.o



clean:
	rm isomorphic_surface database fluid_compression