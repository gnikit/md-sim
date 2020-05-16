SHELL = /bin/bash
include ../Makefile.variable

LIB := ../lib
BIN := ../bin
INCLUDE := ../include

EXCLUDE_FILES := mdmain.cpp main.cpp
FILES := $(filter-out $(EXCLUDE_FILES), $(wildcard *.cpp))
OBJS  := $(patsubst %.cpp, %.o, $(FILES))

all: md

# Adds an xml/rng interface to the MD code
md: libmd.a main.o mdmain.o
	$(CXX) $(CXXFLAGS) mdmain.o main.o $(LINKER) -o $(BIN)/md

# All the MD code stripped from the SPUD interface
libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(LINKER) $< -o $@
	cp $@ $(INCLUDE)/$@


clean:
	$(RM) *.o $(BIN)/md $(LIB)/libmd.a
	$(RM) *.gcda *.gcno *.gcov

