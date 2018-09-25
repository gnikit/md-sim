SHELL = /bin/bash
include ../Makefile.variables
BINDIR = ../bin

# If header files are (and cpp files) are ever separated from
# the libmd.a then -I /dir/to/headers should be used
all: libmd.a

libmd.a: MD.o stat_analysis.o isomorph.o
	ar rcvs libmd.a MD.o stat_analysis.o isomorph.o

MD.o: MD.cpp FileLoading.h
	$(CXX) $(CXXFLAGS) MD.cpp

stat_analysis.o: stat_analysis.cpp
	$(CXX) $(CXXFLAGS) stat_analysis.cpp

isomorph.o: isomorph.cpp
	$(CXX) $(CXXFLAGS) isomorph.cpp

clean:
	rm MD.o stat_analysis.o isomorph.o libmd.a
