SHELL = /bin/bash
include Makefile.variables

RM := rm -rf

default: libmd

all: libmd examples

examples: libmd
	@echo "MAKE MD examples"
	@cd examples && $(MAKE)

libmd:
	@mkdir -p include
	@echo "MAKE lib"
	@cd lib && $(MAKE)
	@echo "MAKE MD src"
	@cd src && $(MAKE)

debug: toolkit
	@echo "DEBUG BUILD"
	@echo "MAKE lib"
	@cd lib && $(MAKE)
	@echo "MAKE MD src"
	@cd src && $(MAKE) debug
	@echo "MAKE examples"
	@cd examples && $(MAKE) debug

schemas: toolkit
	# If the user has installed libspud see fluidity project on github
	@spud-preprocess schemas/main_schema.rnc

toolkit:
	@echo "MAKE tools"
	@cd tools && $(MAKE) scripts

python:
	@echo "MAKE python modules"
	@cd tools/md-tools && pip3 install --user --upgrade -e .

test:
	@echo "Running regression test"
	@cd tests; python run_tests.py

examples: all
	# Do not run the database files
	$(RM) examples/examplebin/*database*
	@cd examples/examplebin; for i in ./*; do echo $$i && ./$$i >> $$i.log; done

clean:
	$(RM) *.log *.txt		# TODO: this will cause an issue if/when we change to cmake
	@echo "Cleaning lib"
	@cd lib && $(MAKE) clean
	@echo "Cleaning MD src and bin"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD Examples src/examples"
	@cd examples && $(MAKE) clean
	@echo "Cleaning bin"
	@cd bin && $(RM) *
	@echo "Cleaning include"
	@cd include && $(RM) *.o

clean_keep_data:
	@echo "Cleaning lib"
	@cd lib && $(MAKE) clean
	@echo "Cleaning MD src and bin"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD Examples src/examples"
	@cd examples && $(MAKE) clean_keep_data
