#define compilers
CC=icc
CXX=icpc
# CXXFLAGS=-c -std=c++17 -pthread -xHOST -use-intel-optimized-headers -m64 -O3 -Wall
# CXXFLAGS=-c -pthread -parallel -mkl=parallel -O3 -xHOST -std=c++17 -use-intel-optimized-headers -m64 -prec-sqrt -prec-div
CXXFLAGS=-c -std=c++17 -pthread -O3 -Wall
SRCDIR=../src

all: project

project: Source.o MD.o stat_analysis.o
	$(CXX) Source.o MD.o stat_analysis.o -o project

Source.o: $(SRCDIR)/Source.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/Source.cpp

MD.o: $(SRCDIR)/MD.cpp 
	$(CXX) $(CXXFLAGS) $(SRCDIR)/MD.cpp

stat_analysis.o: $(SRCDIR)/stat_analysis.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/stat_analysis.cpp

clean:
	rm -rf *o project 

