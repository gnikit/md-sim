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
  std::string dir = "";
  unsigned int steps = 5000;
  bool comp = false;
  unsigned int rdf_bins = 500;
  /* particles per axis */
  unsigned int ppa = 10;
  std::string lattice = "SC";
  bool track_particles = false;
  unsigned int rdf_wait = 2000;
};

struct simulation_type {
  std::string simulation_name = "";
  double density = 0.5;
  double density_final = 2.0;
  double density_inc = 0.5;
  double temperature = 0.5;
  double power = 8;
  double a_cst = 0.5;
  std::string pp_type = "BIP";
};

using namespace tinyxml2;

/*
 * I think there is an easier way to do this although it will not include
 * the detailed error checking that is done now.
 * look at answer:
 * https://stackoverflow.com/questions/54768440/reading-in-all-siblings-elements-with-tinyxml
 */

int load_constructor_options(XMLNode* root, constructor_type& constructor);

int load_simulation_options(XMLNode* root, simulation_type& sim, bool compress);

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

  simulation_type sim;
  XMLCheckResult(load_simulation_options(root, sim, constructor.comp));

  /* Check if the compression flag has been enabled */
  if (constructor.comp) {
    phase_transition run(constructor.dir, constructor.steps, constructor.comp,
                         constructor.rdf_bins, constructor.ppa,
                         constructor.lattice, constructor.track_particles,
                         constructor.rdf_wait);

    run.crystallisation(sim.simulation_name, sim.density, sim.density_final,
                        sim.density_inc, sim.temperature, sim.power, sim.a_cst,
                        sim.pp_type);
  } else {
    MD run(constructor.dir, constructor.steps, constructor.comp,
           constructor.rdf_bins, constructor.ppa, constructor.lattice,
           constructor.track_particles, constructor.rdf_wait);

    run.simulation(sim.simulation_name, sim.density, sim.temperature, sim.power,
                   sim.a_cst, sim.pp_type);
  }
}

int load_constructor_options(XMLNode* root, constructor_type& constructor) {
  /* Pointer to the 1st child of root, constructor */
  XMLNode* input = root->FirstChildElement("constructor");

  /* Create a map to store the data from the XML file */
  std::map<std::string, XMLElement*> vars;

  vars["xml_dir"] = input->FirstChildElement("output_dir");
  vars["xml_steps"] = input->FirstChildElement("steps");
  vars["xml_comp"] = input->FirstChildElement("compression");
  vars["xml_rdf_bins"] = input->FirstChildElement("rdf_bins");
  vars["xml_particles_per_axis"] =
      input->FirstChildElement("particles_per_axis");
  vars["xml_lattice"] = input->FirstChildElement("lattice");
  vars["xml_track_particles"] = input->FirstChildElement("track_particles");
  vars["xml_rdf_wait"] = input->FirstChildElement("rdf_equilibrate");

  /* Check if all results exist and are readbale */
  for (auto const& i : vars) {
    try {
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
  /* Parse output directory */
  constructor.dir = vars["xml_dir"]->GetText();
  /* Parse number of steps */
  error = vars["xml_steps"]->QueryUnsignedText(&constructor.steps);
  XMLCheckResult(error);
  /* Parse compression flag */
  error = vars["xml_comp"]->QueryBoolText(&constructor.comp);
  XMLCheckResult(error);
  /* Parse RDF accuracy */
  error = vars["xml_rdf_bins"]->QueryUnsignedText(&constructor.rdf_bins);
  XMLCheckResult(error);
  /* Parse particles per axis */
  error = vars["xml_particles_per_axis"]->QueryUnsignedText(&constructor.ppa);
  XMLCheckResult(error);
  /* Parse form of lattice */
  constructor.lattice = vars["xml_lattice"]->GetText();
  /* Parse particle tracking flag */
  error =
      vars["xml_track_particles"]->QueryBoolText(&constructor.track_particles);
  XMLCheckResult(error);
  /* Parse RDF equilibration steps */
  error = vars["xml_rdf_wait"]->QueryUnsignedText(&constructor.rdf_wait);
  XMLCheckResult(error);

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
  std::string sim_type = name->Attribute("name");

  /* See if a compressions have been enabled in the schema */
  if (sim_type == "CompressionRun" && compress == true) {
    vars["xml_rho_final"] = input->FirstChildElement("final_density");
    vars["xml_rho_inc"] = input->FirstChildElement("density_increment");
  }

  else if (sim_type == "CompressionRun" && compress == false) {
    std::cerr << "Error: Set compression to true in the schema." << std::endl;
    return -1;
  }

  else if (sim_type != "CompressionRun" && compress == false) {
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

  if (sim_type == "CompressionRun" && compress == true) {
    /* Parse final density for compression */
    error = vars["xml_rho_final"]->QueryDoubleText(&sim.density_final);
    XMLCheckResult(error);
    /* Parse density increment for compression */
    error = vars["xml_rho_inc"]->QueryDoubleText(&sim.density_inc);
    XMLCheckResult(error);
  }

  return 0;
}