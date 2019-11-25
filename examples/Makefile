SHELL = /bin/bash

include ../Makefile.variables
LIB := ../lib
BIN := ../bin
INCLUDE := ../include

all:
	@mkdir -p examplebin example_data
	@for d in $(shell ls -d */); do $(MAKE) -C $${d}; done

debug:
	@mkdir -p examplebin example_data
	@for d in $(shell ls -d */); do $(MAKE) debug -C $${d}; done	

clean:
	@for d in $(shell ls -d */); do $(MAKE) clean -C $${d}; done
	@rm -rf example_data examplebin

clean_keep_data: $(CWD)/*
	@for d in $(shell ls -d */); do $(MAKE) clean_keep_data -C $${d}; done
