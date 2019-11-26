#include <stdlib.h>
#include <tinyxml2.h>
#include <string>
#include <thread>
#include "MD.h"

#ifndef XMLCheckResult
#define XMLCheckResult(a_eResult)     \
  if (a_eResult != XML_SUCCESS) {     \
    printf("Error: %i\n", a_eResult); \
    return a_eResult;                 \
  }
#endif

using namespace tinyxml2;

int main(int argc, char const* argv[]) {
  /* Path and name of the XML schema file */

  XMLDocument schema;
  XMLError eReuslt = schema.LoadFile(argv[1]);
  /* Check if file exists */
  XMLCheckResult(eReuslt);
  /* Pointer to the root of the file */
  XMLNode* root = schema.FirstChildElement();
  /* Check xml file has a first node */
  if (root == nullptr) return XML_ERROR_FILE_READ_ERROR;

  /* Pointer to the 1st child of root, constructor */
  XMLElement* xml_dir =
      root->FirstChildElement("constructor")->FirstChildElement("output_dir");
  XMLElement* xml_steps =
      root->FirstChildElement("constructor")->FirstChildElement("steps");
  XMLElement* xml_comp =
      root->FirstChildElement("constructor")->FirstChildElement("compression");
  XMLElement* xml_rdf_bins =
      root->FirstChildElement("constructor")->FirstChildElement("rdf_bins");
  XMLElement* xml_particles_per_axis =
      root->FirstChildElement("constructor")
          ->FirstChildElement("particles_per_axis");
  XMLElement* xml_lattice =
      root->FirstChildElement("constructor")->FirstChildElement("lattice");
  XMLElement* xml_track_particles = root->FirstChildElement("constructor")
                                        ->FirstChildElement("track_particles");
  XMLElement* xml_rdf_wait = root->FirstChildElement("constructor")
                                 ->FirstChildElement("rdf_equilibrate");

  /* Pointer to the 2nd child of root, simulation_input */
  XMLElement* xml_rho =
      root->FirstChildElement("simulation_input")->FirstChildElement("rho");
  XMLElement* xml_T =
      root->FirstChildElement("simulation_input")->FirstChildElement("T");
  XMLElement* xml_n =
      root->FirstChildElement("simulation_input")->FirstChildElement("n");
  XMLElement* xml_A =
      root->FirstChildElement("simulation_input")->FirstChildElement("A");
  XMLElement* xml_pp =
      root->FirstChildElement("simulation_input")->FirstChildElement("pp");

  /* Check if results exist and are readbale */
  if (xml_dir == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_steps == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_comp == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_rdf_bins == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_particles_per_axis == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_lattice == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_track_particles == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_rdf_wait == nullptr) return XML_ERROR_PARSING_ELEMENT;

  if (xml_rho == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_T == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_n == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_A == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_pp == nullptr) return XML_ERROR_PARSING_ELEMENT;

  /* Convert xml input to usable variables */
  std::string dir, simulation_name, lattice, pair_potential;
  bool comp, track_particles;
  unsigned int steps, n, rdf_bins, particles_per_axis, rdf_wait;
  double rho, T, A;

  dir = xml_dir->GetText();                        // Parse output directory
  eReuslt = xml_steps->QueryUnsignedText(&steps);  // Parse number of steps
  XMLCheckResult(eReuslt);
  eReuslt = xml_comp->QueryBoolText(&comp);  // Parse compression flag
  XMLCheckResult(eReuslt);
  eReuslt = xml_rdf_bins->QueryUnsignedText(&rdf_bins);  // Parse RDF accuracy
  XMLCheckResult(eReuslt);
  eReuslt = xml_particles_per_axis->QueryUnsignedText(
      &particles_per_axis);  // Parse particles per axis
  XMLCheckResult(eReuslt);
  lattice = xml_lattice->GetText();
  eReuslt = xml_track_particles->QueryBoolText(&track_particles);
  XMLCheckResult(eReuslt);
  eReuslt = xml_rdf_wait->QueryUnsignedText(
      &rdf_wait);  // Parse RDF equilibration steps
  XMLCheckResult(eReuslt);

  /* Simulation input */
  simulation_name = xml_dir->GetText();      // Parse output directory
  eReuslt = xml_rho->QueryDoubleText(&rho);  // Parse density
  XMLCheckResult(eReuslt);
  eReuslt = xml_T->QueryDoubleText(&T);  // Parse temperature
  XMLCheckResult(eReuslt);
  eReuslt = xml_n->QueryUnsignedText(&n);  // Parse potential power
  XMLCheckResult(eReuslt);
  eReuslt = xml_A->QueryDoubleText(&A);  // Parse softening parameter
  pair_potential = xml_pp->GetText();    // Parse type of pair potential

  /* Run MD simulation */
  MD run(dir, steps, comp, rdf_bins, particles_per_axis, lattice,
         track_particles, rdf_wait);
  run.simulation(simulation_name, rho, T, n, A, pair_potential);
}
