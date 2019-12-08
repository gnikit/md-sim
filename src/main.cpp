#include <stdlib.h>
#include <tinyxml2.h>
#include <map>
#include <string>
#include <thread>
#include "MD.h"
#include "phase_transition.h"

#ifndef XMLCheckResult
#define XMLCheckResult(a_eResult)     \
  if (a_eResult != XML_SUCCESS) {     \
    printf("Error: %i\n", a_eResult); \
    return a_eResult;                 \
  }
#endif

struct constructor_type {
  unsigned int steps = 5000;
  /* particles per axis */
  unsigned int ppa = 10;
  std::string lattice = "SC";
};

struct options_type {
  std::string dir = ".";
  bool comp = false;
  bool track_particles = false;
  unsigned int rdf_bins = 500;
  unsigned int rdf_wait = 0;
};

struct simulation_type {
  std::string simulation_name = "";
  std::string simulation_type = "";
  double density = 0.5;
  double density_final = 2.0;
  double density_inc = 0.5;
  double temperature = 0.5;
  double power = 8;
  double a_cst = 0.5;
  std::string pp_type = "BIP";
  /* Whether or not to perform a reverse compression run */
  bool reverse_comp = false;
};

struct test_type {
  bool is_testing = false;
};

using namespace tinyxml2;

/*
 * I think there is an easier way to do this although it will not include
 * the detailed error checking that is done now.
 * look at answer:
 * https://stackoverflow.com/questions/54768440/reading-in-all-siblings-elements-with-tinyxml
 */

int get_switch(options_type& options, simulation_type& sim);

int load_constructor_options(XMLNode* root, constructor_type& constructor);

int load_options(XMLNode* root, options_type& options);

int load_simulation_options(XMLNode* root, simulation_type& sim, bool compress);

int load_test_options(XMLNode* root, test_type& test_options);

int main(int argc, char const* argv[]) {
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
  sim.simulation_type = name->Attribute("name");

  /* See if a compressions have been enabled in the schema */
  if (sim.simulation_type == "CompressionRun" && compress == true) {
    vars["xml_rho_final"] = input->FirstChildElement("final_density");
    vars["xml_rho_inc"] = input->FirstChildElement("density_increment");
  }

  else if (sim.simulation_type == "ReverseCompressionRun" && compress == true) {
    vars["xml_rho_final"] = input->FirstChildElement("final_density");
    vars["xml_rho_inc"] = input->FirstChildElement("density_increment");
    sim.reverse_comp = true;
  }

  else if ((sim.simulation_type == "CompressionRun" ||
            sim.simulation_type == "ReverseCompressionRun") &&
           compress == false) {
    std::cerr << "Error: Set constructor/compression to true in the schema."
              << std::endl;
    return -1;
  }

  else if ((sim.simulation_type != "CompressionRun" ||
            sim.simulation_type != "ReverseCompressionRun") &&
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

  if ((sim.simulation_type == "CompressionRun" ||
       sim.simulation_type == "ReverseCompressionRun") &&
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