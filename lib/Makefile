SHELL = /bin/bash
include ../Makefile.variables
spud := ../spud


all:
	@echo "MAKE spud"
	@cd $(spud) && $(MAKE)
	@cp $(spud)/libspud.a .
	@cp $(spud)/include/* ../include
	

clean:
	@echo "Cleaning lib"
	$(RM) *.a
	$(RM) *.o
