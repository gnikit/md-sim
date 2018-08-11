SHELL = /bin/bash
include ../Makefile.variables
LIB =../lib
BIN =../bin

all: project $(LIB)/libmd.a

project: Source.o $(LIB)/libmd.a
		$(CXX) $(FLAGS) Source.o $(LINKER) -o $(BIN)/project

Source.o: Source.cpp
	$(CXX) $(CXXFLAGS) -I$(LIB) Source.cpp -o Source.o

clean:
	rm $(BIN)/project Source.o

