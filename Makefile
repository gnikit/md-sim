SHELL = /bin/bash
include Makefile.variables

RM := rm -rf
SPUD_DIR := spud
CURRENT_DIR := $(shell pwd)

default: libmd schema

all: libmd examples schema


libmd: libspud
	@echo "MAKE MD src"
	@cd src && $(MAKE)

libspud:
	@echo "MAKE libspud"
ifeq ("$(wildcard $(SPUD_DIR)/libspud.a)","")
	@echo "Configuring and compiling libspud"
	@cd $(SPUD_DIR) && ./configure --prefix= --disable-shared
	$(MAKE) -C $(SPUD_DIR) -j
	$(MAKE) -C $(SPUD_DIR) -j install-libspud DESTDIR=..
	$(MAKE) -C $(SPUD_DIR) -j install-spudtools DESTDIR=..
	$(MAKE) -C $(SPUD_DIR) -j install-diamond DESTDIR=../..
	$(MAKE) -C $(SPUD_DIR) -j install-dxdiff DESTDIR=../..
endif
	$(MAKE) -C $(SPUD_DIR) -j
	$(MAKE) -C $(SPUD_DIR) -j install-libspud DESTDIR=..
	$(MAKE) -C $(SPUD_DIR) -j install-spudtools DESTDIR=..
	$(MAKE) -C $(SPUD_DIR) -j install-diamond DESTDIR=../..
	$(MAKE) -C $(SPUD_DIR) -j install-dxdiff DESTDIR=../..

schema:
	# This is a bug fix where because spud-preprocess does not look in the
	# right place for the spud_base.rnc. It ignores the prefix
	sed -i "s+cp /share+cp $(CURRENT_DIR)/share+g" ./bin/spud-preprocess
	./bin/spud-preprocess ./schemas/main_schema.rnc

debug_libmd:
	@echo "DEBUG BUILD"
	@echo "MAKE MD src"
	@cd src && $(MAKE) debug

debug:
	@echo "DEBUG BUILD"
	@echo "MAKE MD src"
	@cd src && $(MAKE) debug
	@echo "MAKE examples"
	@cd examples && $(MAKE) debug

examples: libmd
	@echo "MAKE MD examples"
	@cd examples && $(MAKE)

toolkit:
	@echo "MAKE tools"
	@cd tools && $(MAKE) scripts

python:
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
	$(RM) *.log *.csv
	@echo "Cleaning lib"
	$(RM) lib
	@echo "Cleaning MD src"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD Examples src/examples"
	@cd examples && $(MAKE) clean
	@echo "Cleaning bin"
	$(RM) bin
	@echo "Cleaning share"
	$(RM) share
	@echo "Cleaning include"
	@cd include && $(RM) *.o *.mod spud

distclean: clean
	$(MAKE) -C spud distclean

clean_keep_data:
	@echo "Cleaning lib"
	@cd lib && $(MAKE) clean
	@echo "Cleaning MD src and bin"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD Examples src/examples"
	@cd examples && $(MAKE) clean_keep_data
