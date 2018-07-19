SHELL = /bin/bash
CC  = icc
CXX = g++
# CXXFLAGS=-c -std=c++17 -pthread -xHOST -use-intel-optimized-headers -m64 -O3 -Wall -g
# CXXFLAGS=-c -pthread -parallel -mkl=parallel -O3 -xHOST -std=c++17 -use-intel-optimized-headers -m64 -prec-sqrt -prec-div
CXXFLAGS =-c -std=c++11 -pthread -O3 -Wall
LIBDIR=../../lib
BINDIR=../../bin

all: database iso

database: make_database.o ../MD.o ../stat_analysis.o
	$(CXX) -std=c++11 make_database.o ../MD.o ../stat_analysis.o -lpthread -o database

iso: isomorph_line.o ../MD.o ../isomorph.o
	$(CXX) -std=c++11 -lpthread isomorph_line.o ../MD.o ../isomorph.o -o iso

isomorph_line.o: isomorph_line.cpp
	$(CXX) $(CXXFLAGS) isomorph_line.cpp


make_database.o: make_database.cpp
	$(CXX) $(CXXFLAGS) make_database.cpp


clean:
	rm *.o iso database