SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
SRC = ../src
INC = -I ../lib -I ../src -I ../tools/FileIO
LN = -L ../lib -L ../src -L ../tools/FileIO

LINKER = $(INC) $(LN) -lmd

all: database \
			isomorph_line \
			fluid_compression \
			simple_run \
			visualise_fluid \
			$(LIB)/libmd.a
			

debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: database \
				isomorph_line \
				fluid_compression \
				simple_run \
				visualise_fluid \
				$(LIB)/libmd.a

# TODO: -lstdc++fs has to be added at the end
# TODO: https://stackoverflow.com/questions/33149878/experimentalfilesystem-linker-error
database: database.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) database.o $(LINKER) -o database -lstdc++fs
	
isomorph_line: isomorph_line.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) isomorph_line.o $(LINKER) -o isomorph_line -lstdc++fs

fluid_compression: fluid_compression.o $(LIB)/libmd.a 
	$(CXX) $(FLAGS) fluid_compression.o $(LINKER) -o fluid_compression -lstdc++fs

simple_run: simple_run.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) simple_run.o $(LINKER) -o simple_run -lstdc++fs

visualise_fluid: visualise_fluid.o $(LIB)/libmd.a
	$(CXX) $(FLAGS) visualise_fluid.o $(LINKER) -o visualise_fluid -lstdc++fs


database.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) make_database.cpp -o database.o

isomorph_line.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) isomorph_line.cpp -o isomorph_line.o

fluid_compression.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) fluid_compression.cpp -o fluid_compression.o

simple_run.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) simple_run.cpp -o simple_run.o

visualise_fluid.o:
	$(CXX) $(CXXFLAGS) -I$(SRC) visualise_fluid.cpp -o visualise_fluid.o


clean:
	$(RM) isomorph_line database fluid_compression 
	$(RM) simple_run visualise_fluid *.o 
	$(RM) ./example_data/*.txt
	$(RM) ../data/initial_velocities*
	$(RM) ../data/vx* ../data/vy* ../data/vz*


clean_keep_data:
	$(RM) isomorph_line database fluid_compression 
	$(RM) simple_run visualise_fluid *.o
