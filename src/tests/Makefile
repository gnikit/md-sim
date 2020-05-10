SHELL = /bin/bash
include ../../Makefile.variable

BIN := ../../bin
LIB := ../../lib
INCLUDE := ../../include

INC := -I$(INCLUDE) -I../
LN  := -L$(LIB)
LINKER := $(INC) $(LN) -lmd -lspud -lstdc++fs

FILES := $(wildcard *.cpp)
OBJS  := $(patsubst %.cpp, %.o, $(FILES))

default: test

all: test

test: clean tests-main.o
	$(CXX) $(CXXFLAGS) tests-main.o $(LINKER) -o tests-main
	cp tests-main $(BIN)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(LINKER) $< -o $@

.PHONY: test

clean:
	$(RM) $(OBJS) tests-main
