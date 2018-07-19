#define compilers
CC=icc
CXX=g++
# CXXFLAGS=-c -std=c++17 -pthread -xHOST \
# 				 -use-intel-optimized-headers -m64 -O3 -Wall
# CXXFLAGS=-c -pthread -parallel -mkl=parallel -O3 -xHOST -std=c++17 -use-intel-optimized-headers -m64 -prec-sqrt -prec-div
CXXFLAGS=-c -std=c++11 -pthread -O0 -Wall
SRCDIR=../src
LIBDIR=../lib
objects = Source.o MD.o stat_analysis.o \
					isomorph.o
all: project

project: $(objects)
	$(CXX) -std=c++11 $(objects) -lpthread -o project

Source.o: $(SRCDIR)/Source.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/Source.cpp

MD.o: $(SRCDIR)/MD.h $(LIBDIR)/FileLoading.h
	$(CXX) $(CXXFLAGS) $(SRCDIR)/MD.cpp

isomorph.o: $(SRCDIR)/isomorph.h
	$(CXX) $(CXXFLAGS) $(SRCDIR)/isomorph.cpp

stat_analysis.o: $(SRCDIR)/stat_analysis.h
	$(CXX) $(CXXFLAGS) $(SRCDIR)/stat_analysis.cpp

clean:
	rm $(objects) project 

