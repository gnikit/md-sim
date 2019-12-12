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
  std::string simulation_form = "";
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

int mdmain(int argc, char const* argv[]);

int get_switch(options_type& options, simulation_type& sim);

int load_constructor_options(XMLNode* root, constructor_type& constructor);

int load_options(XMLNode* root, options_type& options);

int load_simulation_options(XMLNode* root, simulation_type& sim, bool compress);

int load_test_options(XMLNode* root, test_type& test_options);
