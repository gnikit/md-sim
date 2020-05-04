SHELL = /bin/bash
include ../Makefile.variables
LIB := ../lib
BIN := ../bin
INCLUDE := ../include

FILES := $(wildcard *.cpp)
OBJS  := $(patsubst %.cpp, %.o, $(FILES))

all: md libmd.a

debug: CXXFLAGS = $(DCXXFLAGS)
debug: md libmd.a

md: libmd.a main.o 
	$(CXX) $(CXXFLAGS) main.o $(LINKER) -o $(BIN)/md

libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(LINKER) $< -o $@
	cp $@ $(INCLUDE)/$@


clean:
	$(RM) $(OBJS) $(BIN)/md $(LIB)/libmd.a

