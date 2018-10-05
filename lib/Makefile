SHELL = /bin/bash
include ../Makefile.variables
tinyxml = ../tools/tinyxml2


all:
	@echo "MAKE tools/TinyXML2"
	@cd $(tinyxml) && $(MAKE)
	@mv $(tinyxml)/libtinyxml2.a .
	@mv $(tinyxml)/xmltest ../bin
	

clean:
	$(RM) libmd.a libtinyxml2.a
	$(RM) *.o
