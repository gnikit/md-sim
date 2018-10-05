SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
OBJS = MD.o stat_analysis.o isomorph.o


all: libmd.a md

md: Source.o libmd.a
		$(CXX) $(FLAGS) Source.o $(LINKER) -o $(BIN)/md

libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

MD.o: MD.cpp
	$(CXX) $(CXXFLAGS) $(LINKER) MD.cpp

stat_analysis.o: stat_analysis.cpp
	$(CXX) $(CXXFLAGS) stat_analysis.cpp

isomorph.o: isomorph.cpp
	$(CXX) $(CXXFLAGS) isomorph.cpp

Source.o: Source.cpp 
	$(CXX) $(CXXFLAGS) $(INC) Source.cpp -o Source.o


clean:
	$(RM) *.o $(BIN)/md $(LIB)/libmd.a

