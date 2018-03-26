#define compiler
CC=gcc
CXX=g++
CXXFLAGS=-c -pthread -parallel -mkl=parallel -O3 -xHOST -std=c++17 -use-intel-optimized-headers -m64 -prec-sqrt -prec-div
SRCDIR=../src

all: project

project: main.o source_files.o
	$(CXX) main.o source_files.o -o project
	# not sure if -o project needs to be the same as project:

main.o: $(SRCDIR)/Source.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/Source.cpp

source_files.o: $(SRCDIR)/MD.cpp $(SRCDIR)/stat_analysis.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/MD.cpp $(SRCDIR)/stat_analysis.cpp -o source_files.o

clean:
	rm -rf *o project 

