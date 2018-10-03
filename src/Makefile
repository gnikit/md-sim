SHELL = /bin/bash
include ../Makefile.variables
LIB = ../lib
BIN = ../bin
# LN = -I../FileIO -L../FileIO
OBJS = MD.o stat_analysis.o isomorph.o

# If header files are (and cpp files) are ever separated from
# the libmd.a then -I /dir/to/headers should be used
all: libmd.a md

md: Source.o libmd.a
		$(CXX) $(FLAGS) Source.o $(LINKER) -o $(BIN)/md

libmd.a: $(OBJS)
	ar rcvs $(LIB)/libmd.a $(OBJS)

MD.o: MD.cpp
	$(CXX) $(CXXFLAGS) $(LINKER) MD.cpp

stat_analysis.o: stat_analysis.cpp
	$(CXX) $(CXXFLAGS) stat_analysis.cpp

isomorph.o: isomorph.cpp
	$(CXX) $(CXXFLAGS) isomorph.cpp

Source.o: Source.cpp
	$(CXX) $(CXXFLAGS) -I$(LIB) Source.cpp -o Source.o

clean:
	$(RM) *.o $(BIN)/md $(LIB)/libmd.a

