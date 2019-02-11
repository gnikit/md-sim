SHELL = /bin/bash

# TODO: -lstdc++fs has to be added at the end
# TODO: https://stackoverflow.com/questions/33149878/experimentalfilesystem-linker-error

all:
	@echo "3D_visualisation"
	@cd 3D_visualisation && $(MAKE)
	@echo "BIP_database"
	@cd BIP_database && $(MAKE)
	@echo "BIP_isomorphic_surfaces"
	@cd BIP_isomorphic_surfaces && $(MAKE)
	@echo "BIP_phase_diagram"
	@cd BIP_phase_diagram && $(MAKE)
	@echo "BIP_single"
	@cd BIP_single && $(MAKE)
	@echo "EXP_database"
	@cd EXP_database && $(MAKE)


clean:
	@echo "3D_visualisation"
	@cd 3D_visualisation && $(MAKE) clean
	@echo "BIP_database"
	@cd BIP_database && $(MAKE) clean
	@echo "BIP_isomorphic_surfaces"
	@cd BIP_isomorphic_surfaces && $(MAKE) clean
	@echo "BIP_phase_diagram"
	@cd BIP_phase_diagram && $(MAKE) clean
	@echo "BIP_single"
	@cd BIP_single && $(MAKE) clean
	@echo "EXP_database"
	@cd EXP_database && $(MAKE) clean

clean_keep_data:
	@echo "3D_visualisation"
	@cd 3D_visualisation && $(MAKE) clean_keep_data
	@echo "BIP_database"
	@cd BIP_database && $(MAKE) clean_keep_data
	@echo "BIP_isomorphic_surfaces"
	@cd BIP_isomorphic_surfaces && $(MAKE) clean_keep_data
	@echo "BIP_phase_diagram"
	@cd BIP_phase_diagram && $(MAKE) clean_keep_data
	@echo "BIP_single"
	@cd BIP_single && $(MAKE) clean_keep_data
	@echo "EXP_database"
	@cd EXP_database && $(MAKE) clean_keep_data
