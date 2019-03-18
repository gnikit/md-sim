SHELL = /bin/bash
include ../Makefile.variables
tinyxml := ../tools/tinyxml2


all:
	@echo "MAKE tools/TinyXML2"
	@cd $(tinyxml) && $(MAKE) staticlib
	@mv $(tinyxml)/libtinyxml2.a .
	@mv $(tinyxml)/tinyxml2.o ../include
	

clean:
	@echo "Cleaning TinyXML2"
	@cd $(tinyxml) && $(MAKE) clean
	@echo "Cleaning lib"
	$(RM) *.a
	$(RM) *.o
