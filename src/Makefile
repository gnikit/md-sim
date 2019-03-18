SHELL = /bin/bash
include ../Makefile.variables
LIB := ../lib
BIN := ../bin
INCLUDE := ../include

# FILES := $(filter-out Source.cpp , $(basename $(wildcard *.cpp)))
FILES := $(basename $(wildcard *.cpp))
OBJS = MD.o isomorph.o md_pair_potentials.o stat_file_logger.o phase_transition.o

all: md libmd.a

debug: CXXFLAGS = $(DCXXFLAGS)
debug: md libmd.a

md: libmd.a Source.o 
	$(CXX) $(CXXFLAGS) Source.o $(LINKER) -o $(BIN)/md

libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(LINKER) $<
	cp $@ $(INCLUDE)/$@


clean:
	$(RM) *.o $(BIN)/md $(LIB)/libmd.a

