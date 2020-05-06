SHELL = /bin/bash
include Makefile.variable

RM := rm -rf
SPUD_DIR := spud
CURRENT_DIR := $(shell pwd)

default: libmd schemas

all: libmd examples schemas

libmd: libspud fileio
	@echo "MAKE MD src"
	$(MAKE) -C src

libspud:
ifeq ("$(wildcard $(SPUD_DIR)/libspud.a)","")
	@echo "Configuring libspud"
	@cd $(SPUD_DIR) && ./configure --prefix= --disable-shared
endif
	@echo "MAKE libspud"
	$(MAKE) -C $(SPUD_DIR)
	$(MAKE) -C $(SPUD_DIR) install-libspud DESTDIR=..
	$(MAKE) -C $(SPUD_DIR) install-spudtools DESTDIR=..
	$(MAKE) -C $(SPUD_DIR) install-diamond DESTDIR=../..
	$(MAKE) -C $(SPUD_DIR) install-dxdiff DESTDIR=../..

fileio:
	cp tools/FileIO/FileIO.h include/

schemas: libspud
	# This is a bug fix where because spud-preprocess does not look in the
	# right place for the spud_base.rnc. It ignores the prefix
	sed -i "s+cp /share+cp $(CURRENT_DIR)/share+g" ./bin/spud-preprocess
	./bin/spud-preprocess ./schemas/main_schema.rnc

.PHONY: scehmas

examples: libmd
	@echo "MAKE MD examples"
	$(MAKE) -C examples

toolkit:
	@echo "MAKE tools"
	$(MAKE) -C tools scripts

python-md-tools:
	@echo "MAKE python modules"
	@cd tools/md-tools && pip3 install --user --upgrade -e .

test: libmd
	@echo "Running regression test"
	@cd tests; python3 run_tests.py

test_examples: libmd
	# Do not run the database files
	$(RM) examples/examplebin/*database*
	@cd examples/examplebin; for i in ./*; do echo $$i && ./$$i >> $$i.log; done

clean:
	@echo "Cleaning lib"
	$(RM) lib
	@echo "Cleaning bin"
	$(RM) bin
	@echo "Cleaning share"
	$(RM) share
	@echo "Cleaning MD src"
	$(MAKE) -C src clean
	@echo "Cleaning examples"
	$(MAKE) -C examples clean
	@echo "Cleaning include"
	@cd include && $(RM) *.o *.mod spud

distclean: clean
	$(MAKE) -C spud distclean
