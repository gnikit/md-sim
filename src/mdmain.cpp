#include "mdmain.h"

int mdmain(int argc, char const* argv[]) {
  /* Path and name of the XML schema file */
  XMLDocument schema;
  XMLError error = schema.LoadFile(argv[1]);
  /* Check if file exists */
  XMLCheckResult(error);

  /* Pointer to the root of the file */
  XMLNode* root = schema.FirstChildElement();
  /* Check xml file has a first node */
  if (root == nullptr) return XML_ERROR_FILE_READ_ERROR;

  constructor_type constructor;
  XMLCheckResult(load_constructor_options(root, constructor));

  options_type options;
  XMLCheckResult(load_options(root, options));

  simulation_type sim;
  XMLCheckResult(load_simulation_options(root, sim, options.comp));

  test_type test;
  XMLCheckResult(load_test_options(root, test));

  /* Run a version of the simulation */
  switch (get_switch(options, sim)) {
    case 0: {
      MD run(constructor.steps, constructor.ppa, constructor.lattice);

      run.load_options(options.dir, options.comp, options.track_particles,
                       options.rdf_bins, options.rdf_wait);

      /* If we are testing this changes to true and fixes the random seed */
      run.enable_testing(test.is_testing);

      run.simulation(sim.simulation_name, sim.density, sim.temperature,
                     sim.power, sim.a_cst, sim.pp_type);
      break;
    }

    case 1: {
      phase_transition run(constructor.steps, constructor.ppa,
                           constructor.lattice);

      run.load_options(options.dir, options.comp, options.track_particles,
                       options.rdf_bins, options.rdf_wait);

      /* If we are testing this changes to true and fixes the random seed */
      run.enable_testing(test.is_testing);

      run.crystallisation(sim.simulation_name, sim.density, sim.density_final,
                          sim.density_inc, sim.temperature, sim.power,
                          sim.a_cst, sim.pp_type);
      break;
    }

    case 2: {
      phase_transition run(constructor.steps, constructor.ppa,
                           constructor.lattice);

      run.load_options(options.dir, options.comp, options.track_particles,
                       options.rdf_bins, options.rdf_wait);

      /* If we are testing this changes to true and fixes the random seed */
      run.enable_testing(test.is_testing);

      run.run_backwards(sim.simulation_name, sim.density, sim.density_final,
                        sim.density_inc, sim.temperature, sim.power, sim.a_cst,
                        sim.pp_type);
      break;
    }

    default:
      break;
  }
  return 0;
}

int get_switch(options_type& options, simulation_type& sim) {
  /*
   * Returns:
   * 0: NormalRun
   * 1: CompressionRun
   * 2: ReverseCompressionRun
   */
  return options.comp + sim.reverse_comp;
}

int load_constructor_options(XMLNode* root, constructor_type& constructor) {
  /* Pointer to the 1st child of root, constructor */
  XMLNode* input = root->FirstChildElement("constructor");

  /* Create a map to store the data from the XML file */
  std::map<std::string, XMLElement*> vars;

  vars["xml_steps"] = input->FirstChildElement("steps");
  vars["xml_particles_per_axis"] =
      input->FirstChildElement("particles_per_axis");
  vars["xml_lattice"] = input->FirstChildElement("lattice");

  /* Check if all the mandatory results exist and are readbale */
  for (auto const& i : vars) {
    try {
      if (i.second == nullptr && i.first.find("xml_opt_") == std::string::npos)
        throw "XML node " + i.first + " missing from the input file.";

    } catch (const char* msg) {
      std::cerr << "Error: " << msg << std::endl;
      exit(1);
    }
  }

  /* Cast XML input to data types with error checking */
  /* The error variable will return an integer with the XML node erroring */
  XMLError error;
  /* Parse number of steps */
  error = vars["xml_steps"]->QueryUnsignedText(&constructor.steps);
  XMLCheckResult(error);
  /* Parse particles per axis */
  error = vars["xml_particles_per_axis"]->QueryUnsignedText(&constructor.ppa);
  XMLCheckResult(error);
  /* Parse form of lattice */
  constructor.lattice = vars["xml_lattice"]->GetText();

  return 0;
}

int load_options(XMLNode* root, options_type& options) {
  /* Pointer to the 1st child of root, options */
  // todo: check if options are present
  XMLNode* input = root->FirstChildElement("options");

  /* Create a map to store the data from the XML file */
  std::map<std::string, XMLElement*> vars;

  /* Optional parameters start with xml_opt_ when parsed */
  vars["xml_opt_dir"] = input->FirstChildElement("output_dir");
  vars["xml_opt_comp"] = input->FirstChildElement("compression");
  vars["xml_opt_track_particles"] = input->FirstChildElement("track_particles");
  vars["xml_opt_rdf_bins"] = input->FirstChildElement("rdf_bins");
  vars["xml_opt_rdf_wait"] = input->FirstChildElement("rdf_equilibrate");

  /* Cast XML input to data types with error checking */
  /* The error variable will return an integer with the XML node erroring */
  XMLError error;

  /* Do not parse the empty optional arguments */
  /* Parse output directory */
  if (vars["xml_opt_dir"] != nullptr) {
    options.dir = vars["xml_opt_dir"]->GetText();
  }
  /* Parse compression flag */
  if (vars["xml_opt_comp"] != nullptr) {
    error = vars["xml_opt_comp"]->QueryBoolText(&options.comp);
    XMLCheckResult(error);
  }
  /* Parse particle tracking flag */
  if (vars["xml_opt_track_particles"] != nullptr) {
    error = vars["xml_opt_track_particles"]->QueryBoolText(
        &options.track_particles);
    XMLCheckResult(error);
  }
  /* Parse RDF accuracy */
  if (vars["xml_opt_rdf_bins"] != nullptr) {
    error = vars["xml_opt_rdf_bins"]->QueryUnsignedText(&options.rdf_bins);
    XMLCheckResult(error);
  }
  /* Parse RDF equilibration steps */
  if (vars["xml_opt_rdf_wait"] != nullptr) {
    error = vars["xml_opt_rdf_wait"]->QueryUnsignedText(&options.rdf_wait);
    XMLCheckResult(error);
  }

  return 0;
}

int load_simulation_options(XMLNode* root, simulation_type& sim,
                            bool compress = false) {
  /* Pointer to the 2nd child of root, simulation_input */
  XMLNode* input = root->FirstChildElement("simulation_input");

  /* Create a map to store the data from the XML file */
  std::map<std::string, XMLElement*> vars;

  /* Get the type of simulation */
  XMLElement* name = input->ToElement();
  sim.simulation_form = name->Attribute("name");

  /* See if a compressions have been enabled in the schema */
  if (sim.simulation_form == "CompressionRun" && compress == true) {
    vars["xml_rho_final"] = input->FirstChildElement("final_density");
    vars["xml_rho_inc"] = input->FirstChildElement("density_increment");
  }

  else if (sim.simulation_form == "ReverseCompressionRun" && compress == true) {
    vars["xml_rho_final"] = input->FirstChildElement("final_density");
    vars["xml_rho_inc"] = input->FirstChildElement("density_increment");
    sim.reverse_comp = true;
  }

  else if ((sim.simulation_form == "CompressionRun" ||
            sim.simulation_form == "ReverseCompressionRun") &&
           compress == false) {
    std::cerr << "Error: Set constructor/compression to true in the schema."
              << std::endl;
    return -1;
  }

  else if ((sim.simulation_form != "CompressionRun" ||
            sim.simulation_form != "ReverseCompressionRun") &&
           compress == true) {
    std::cerr << "Error: Set simulation type to CompressionRun in the schema."
              << std::endl;
    return -1;
  }

  vars["xml_simulation_name"] = input->FirstChildElement("simulation_name");
  vars["xml_rho"] = input->FirstChildElement("rho");
  vars["xml_T"] = input->FirstChildElement("T");
  vars["xml_n"] = input->FirstChildElement("n");
  vars["xml_A"] = input->FirstChildElement("A");
  vars["xml_pp"] = input->FirstChildElement("pp");

  /* Check if all results exist and are readbale */
  for (auto const& i : vars) {
    try {
      // todo: check if string in vector of strings then continue
      if (i.second == nullptr)
        throw "XML node " + i.first + " missing from the input file.";

    } catch (const char* msg) {
      std::cerr << "Error: " << msg << std::endl;
      exit(1);
    }
  }

  /* Cast XML input to data types with error checking */
  /* The error variable will return an integer with the XML node erroring */
  XMLError error;
  /* Parse simulation name */
  sim.simulation_name = vars["xml_simulation_name"]->GetText();
  /* Parse density and check for errors */
  error = vars["xml_rho"]->QueryDoubleText(&sim.density);
  XMLCheckResult(error);
  /* Parse temperature and check for errors */
  error = vars["xml_T"]->QueryDoubleText(&sim.temperature);
  XMLCheckResult(error);
  /* Parse potential power and check for errors */
  error = vars["xml_n"]->QueryDoubleText(&sim.power);
  XMLCheckResult(error);
  /* Parse softening parameter and check for errors */
  error = vars["xml_A"]->QueryDoubleText(&sim.a_cst);
  XMLCheckResult(error);
  /* Parse type of pair potential */
  sim.pp_type = vars["xml_pp"]->GetText();

  if ((sim.simulation_form == "CompressionRun" ||
       sim.simulation_form == "ReverseCompressionRun") &&
      compress == true) {
    /* Parse final density for compression */
    error = vars["xml_rho_final"]->QueryDoubleText(&sim.density_final);
    XMLCheckResult(error);
    /* Parse density increment for compression */
    error = vars["xml_rho_inc"]->QueryDoubleText(&sim.density_inc);
    XMLCheckResult(error);
  }

  return 0;
}

int load_test_options(XMLNode* root, test_type& test_options) {
  XMLElement* xml_is_testing = root->FirstChildElement("enable_testing");

  XMLError error;
  /* Parse the testing flag from the schema */
  error = xml_is_testing->QueryBoolText(&test_options.is_testing);
  XMLCheckResult(error);

  return 0;
}