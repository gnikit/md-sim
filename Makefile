SHELL = /bin/bash
include Makefile.variables

RM := rm -rf

all: libmd examples schemas

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
	@bin/generate_xml_file -d . -s 5000 -c false -r 500 -p 10 -l SC -t false -w 2000 -N simple_sc_run_ -R 0.5 -T 0.5 -n 8 -A 0.5 -P BIP -o schemas/bip_sc_run_input >/dev/null
	@bin/generate_xml_file -d . -s 5000 -c false -r 500 -p 10 -l FCC -t false -w 2000 -N simple_fcc_run_ -R 0.5 -T 0.5 -n 8 -A 0.5 -P BIP -o schemas/bip_fcc_run_input >/dev/null
	@bin/generate_xml_file -d . -s 5000 -c false -r 500 -p 10 -l BCC -t false -w 2000 -N simple_bcc_run_ -R 0.5 -T 0.5 -n 8 -A 0.5 -P BIP -o schemas/bip_bcc_run_input >/dev/null
	@bin/generate_xml_file -d . -s 5000 -c false -r 500 -p 10 -l FCC -t false -w 2000 -N simple_bcc_run_ -R 0.5 -T 0.5 -n 8 -A 0.5 -P GCM -o schemas/gcm_fcc_run_input >/dev/null
	@bin/generate_xml_file -d . -s 5000 -c false -r 500 -p 10 -l BCC -t false -w 2000 -N simple_bcc_run_ -R 0.5 -T 0.5 -n 8 -A 0.5 -P GCM -o schemas/gcm_bcc_run_input >/dev/null

toolkit:
	@echo "MAKE tools"
	@cd tools && $(MAKE) scripts

python:
	@echo "MAKE python modules"
	@cd tools/md-tools && pip3 install --user --upgrade -e .

test: all
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
