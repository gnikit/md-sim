SHELL = /bin/bash

CC  = icc
CXX = g++
# CXXFLAGS=-c -std=c++17 -pthread -xHOST -use-intel-optimized-headers -m64 -O3 -Wall -g
# CXXFLAGS=-c -pthread -parallel -mkl=parallel -O3 -xHOST -std=c++17 -use-intel-optimized-headers -m64 -prec-sqrt -prec-div
CXXFLAGS =-c -std=c++11 -pthread -O3 -Wall
LIBDIR=../lib
BINDIR=../bin
objects = Source.o \
					MD.o \
					stat_analysis.o \
					isomorph.o

all: project

project: $(objects)
	$(CXX) -std=c++11 $(objects) -lpthread -o $(BINDIR)/project

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



