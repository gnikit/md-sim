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
  XMLElement* xml_dir = root->FirstChildElement("constructor")->FirstChildElement("output_dir");
  XMLElement* xml_steps = root->FirstChildElement("constructor")->FirstChildElement("steps");
  XMLElement* xml_comp = root->FirstChildElement("constructor")->FirstChildElement("compression");

  /* Pointer to the 2nd child of root, simulation_input */
  XMLElement* xml_rho = root->FirstChildElement("simulation_input")->FirstChildElement("rho");
  XMLElement* xml_T = root->FirstChildElement("simulation_input")->FirstChildElement("T");
  XMLElement* xml_n = root->FirstChildElement("simulation_input")->FirstChildElement("n");
  XMLElement* xml_A = root->FirstChildElement("simulation_input")->FirstChildElement("A");

  /* Check if reasults exist and are readbale */
  if (xml_dir == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_steps == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_rho == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_T == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_n == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_A == nullptr) return XML_ERROR_PARSING_ELEMENT;
  if (xml_comp == nullptr) return XML_ERROR_PARSING_ELEMENT;

  /* Convert xml input to usable variables */
  std::string dir;
  bool comp;
  unsigned int steps, n;
  double rho, T, A;

  dir = xml_dir->GetText();                        // Parse output directory
  eReuslt = xml_steps->QueryUnsignedText(&steps);  // Parse number of steps
  XMLCheckResult(eReuslt);
  eReuslt = xml_comp->QueryBoolText(&comp);  // Parse compression flag
  XMLCheckResult(eReuslt);
  eReuslt = xml_rho->QueryDoubleText(&rho);  // Parse density
  XMLCheckResult(eReuslt);
  eReuslt = xml_T->QueryDoubleText(&T);  // Parse temperature
  XMLCheckResult(eReuslt);
  eReuslt = xml_n->QueryUnsignedText(&n);  // Parse potential power
  XMLCheckResult(eReuslt);
  eReuslt = xml_A->QueryDoubleText(&A);  // Parse softening parameter

  /* Run MD simulation */
  MD run(dir, steps);
  run.Simulation(rho, T, n, A);
}
