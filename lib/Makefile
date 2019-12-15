SHELL = /bin/bash
include ../Makefile.variables
spud := ../spud


all:
	@echo "MAKE spud"
ifeq ("$(wildcard $(spud)/libspud.a)","")
	@echo "Configuring and compiling libspud"
	@cd $(spud) && ./configure && $(MAKE)
endif
	@cp $(spud)/libspud.a .
	@cp $(spud)/include/* ../include
	

clean:
	@echo "Cleaning lib"
	$(RM) *.a
	$(RM) *.o
