SHELL = /bin/bash


all:
	@echo "MAKE MD src"
	@cd src && $(MAKE)
	@echo "MAKE MD examples"
	cd src/examples && $(MAKE)

clean:
	@echo "Cleaning MD ./src"
	@cd  src && $(MAKE) clean
	@echo "Cleaning MD ./bin"
	@cd bin && $(MAKE) clean
	@echo "Cleaning MD Examples ./src/examples"
	@cd src/examples && $(MAKE) clean
