SHELL = /bin/bash


all:
	@echo "MAKE tools/TinyXML2"
	@cd tools/tinyxml2 && $(MAKE) staticlib
	@echo "MAKE lib"
	@cd lib && $(MAKE)
	@echo "MAKE MD src"
	@cd src && $(MAKE)
	@echo "MAKE MD examples"
	@cd examples && $(MAKE)
	@echo "MAKE tools"
	@cd tools && $(MAKE)

clean:
	@echo "Cleaning lib"
	@cd lib && $(MAKE) clean
	@echo "Cleaning TinyXML2"
	@cd tools/tinyxml2 && $(MAKE) clean
	@echo "Cleaning MD src and bin"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD Examples src/examples"
	@cd examples && $(MAKE) clean
	@echo "MAKE tools"
	@cd tools && $(MAKE) clean

debug:
	@echo "DEBUG BUILD"
	@echo "MAKE tools/TinyXML2"
	@cd tools/tinyxml2 && $(MAKE) staticlib
	@echo "MAKE lib"
	@cd lib && $(MAKE)
	@echo "MAKE MD src"
	@cd src && $(MAKE) debug
	@echo "MAKE MD examples"
	@cd examples && $(MAKE)
	@echo "MAKE tools"
	@cd tools && $(MAKE)

clean_keep_data:
	@echo "Cleaning lib"
	@cd lib && $(MAKE) clean
	@echo "Cleaning TinyXML2"
	@cd tools/tinyxml2 && $(MAKE) clean
	@echo "Cleaning MD src and bin"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD Examples src/examples"
	@cd examples && $(MAKE) clean_keep_data
	@echo "MAKE tools"
	@cd tools && $(MAKE) clean
