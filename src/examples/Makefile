SHELL = /bin/bash
include ../Makefile.variables
LIBDIR=../../lib
BINDIR=../../bin
# TODO: This will obviously not work if isomorph.o, MD.o stat_analysis.o do not exist
# 			at this point 'make all' has to be called and 'database' and 'iso' cannot be called independently
all: database iso compression

database: make_database.o ../MD.o ../stat_analysis.o
	$(CXX) -std=$(STD) make_database.o ../MD.o ../stat_analysis.o -pthread -o database

iso: isomorph_line.o ../MD.o ../isomorph.o
	$(CXX) -std=$(STD) -pthread isomorph_line.o ../MD.o ../isomorph.o -o iso

compression: fluid_compression.o ../MD.o 
	$(CXX) -std=$(STD) -pthread fluid_compression.o ../MD.o -o compression

isomorph_line.o: 
	$(CXX) $(CXXFLAGS) isomorph_line.cpp

make_database.o: 
	$(CXX) $(CXXFLAGS) make_database.cpp

fluid_compression.o:
	$(CXX) $(CXXFLAGS) fluid_compression.cpp

clean:
	rm *.o iso database