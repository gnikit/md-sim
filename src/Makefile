SHELL = /bin/bash
include Makefile.variables

LIBDIR=../lib
BINDIR=../bin
objects = Source.o \
	  MD.o \
	  stat_analysis.o \
	  isomorph.o

all: project

project: $(objects)
	$(CXX) -std=$(STD) $(objects) -pthread -o $(BINDIR)/project

Source.o: Source.cpp
	$(CXX) $(CXXFLAGS) Source.cpp

MD.o: MD.h $(LIBDIR)/FileLoading.h
	$(CXX) $(CXXFLAGS) MD.cpp

isomorph.o: isomorph.h
	$(CXX) $(CXXFLAGS) isomorph.cpp

stat_analysis.o: stat_analysis.h
	$(CXX) $(CXXFLAGS) stat_analysis.cpp

clean:
	rm $(objects)



