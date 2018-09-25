SHELL = /bin/bash
include ../Makefile.variables
LIB =../lib
BIN =../bin

all: md $(LIB)/libmd.a

md: Source.o $(LIB)/libmd.a
		$(CXX) $(FLAGS) Source.o $(LINKER) -o $(BIN)/md

Source.o: Source.cpp
	$(CXX) $(CXXFLAGS) -I$(LIB) Source.cpp -o Source.o

clean:
	rm $(BIN)/md Source.o

