SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
OBJS = MD.o isomorph.o md_algorithms.o

# The LINKER and the INC variables contain invalid, or unnecassary paths
# like ../src or ../tools/FileIO for the case of Source.o 
#which cause the compilation to be slightly slower. 
# The trade off is global-management of the linker from Makefile.variables

all: md libmd.a

debug: CXXFLAGS = $(DCXXFLAGS)
debug: FLAGS = $(DFLAGS)
debug: md libmd.a

md: Source.o libmd.a
		$(CXX) $(FLAGS) Source.o $(LINKER) -o $(BIN)/md

libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

MD.o: MD.cpp
	$(CXX) $(CXXFLAGS) $(LINKER) MD.cpp

md_algorithms.o: md_algorithms.cpp
	$(CXX) $(CXXFLAGS) $(LINKER) md_algorithms.cpp

isomorph.o: isomorph.cpp
	$(CXX) $(CXXFLAGS) isomorph.cpp

Source.o: Source.cpp 
	$(CXX) $(CXXFLAGS) $(INC) Source.cpp


clean:
	$(RM) *.o $(BIN)/md $(LIB)/libmd.a

