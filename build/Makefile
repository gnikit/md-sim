#define compilers
SHELL = /bin/bash
PACKAGE_NAME = MDsim
CC  = icc
CXX = g++
# CXXFLAGS=-c -std=c++17 -pthread -xHOST \
# 				 -use-intel-optimized-headers -m64 -O3 -Wall
# CXXFLAGS=-c -pthread -parallel -mkl=parallel -O3 -xHOST -std=c++17 -use-intel-optimized-headers -m64 -prec-sqrt -prec-div
CXXFLAGS = -c -std=c++11 -pthread -O3 -Wall
SRCDIR   = ../src
LIBDIR   = ../lib
BINDIR   = ../bin
objects  = MD.o stat_analysis.o isomorph.o
source   = Source.o
data     = build_database.o	

all: project database

project: $(objects) $(source)
	$(CXX) -std=c++11 $(source) $(objects) -lpthread -o $(BINDIR)/project

database: $(objects) $(data)
	$(CXX) -std=c++11 $(data) $(objects) -lpthread -o $(BINDIR)/database

build_database.o: $(SRCDIR)/build_database.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/build_database.cpp

Source.o: $(SRCDIR)/Source.cpp
	$(CXX) $(CXXFLAGS) $(SRCDIR)/Source.cpp

MD.o: $(SRCDIR)/MD.h $(LIBDIR)/FileLoading.h
	$(CXX) $(CXXFLAGS) $(SRCDIR)/MD.cpp

isomorph.o: $(SRCDIR)/isomorph.h
	$(CXX) $(CXXFLAGS) $(SRCDIR)/isomorph.cpp

stat_analysis.o: $(SRCDIR)/stat_analysis.h
	$(CXX) $(CXXFLAGS) $(SRCDIR)/stat_analysis.cpp

clean:
	rm $(objects) $(source) ../bin/project ../bin/database

