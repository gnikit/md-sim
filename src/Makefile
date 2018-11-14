SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
OBJS = MD.o isomorph.o md_algorithms.o


all: libmd.a md

md: Source.o libmd.a
		$(CXX) $(FLAGS) Source.o $(LINKER) -o $(BIN)/md

libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

# TODO: maybe need to include md_algorithms.cpp
MD.o: MD.cpp
	$(CXX) $(CXXFLAGS) $(LINKER) MD.cpp

md_algorithms.o: md_algorithms.cpp
	$(CXX) $(CXXFLAGS) $(LINKER) md_algorithms.cpp

isomorph.o: isomorph.cpp
	$(CXX) $(CXXFLAGS) isomorph.cpp

Source.o: Source.cpp 
	$(CXX) $(CXXFLAGS) $(INC) Source.cpp -o Source.o


clean:
	$(RM) *.o $(BIN)/md $(LIB)/libmd.a

