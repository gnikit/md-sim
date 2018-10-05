SHELL = /bin/bash
include ../Makefile.variables
LIB = ../bin
tinyxml = ../tools/tinyxml2


all:
	@echo "MAKE tools/TinyXML2"
	@cd $(tinyxml) && $(MAKE)
	
	
clean:
	$(RM) libmd.a libtinyxml2.a
	$(RM) *.o
