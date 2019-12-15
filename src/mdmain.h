#include <stdlib.h>
#include <any>
#include <map>
#include <string>
#include <thread>
#include "MD.h"
#include "phase_transition.h"
#include "spud"

struct constructor_type {
  size_t steps = 0;
  std::vector<size_t> particles;
  std::string lattice = "";
  double random_variance = 0;
};

struct options_type {
  std::string dir = ".";
  bool track_particles = false;
  size_t rdf_bins = 0;
  size_t rdf_wait = 0;
};

struct simulation_type {
  std::string simulation_type = "";
  std::string simulation_name = "";
  double density = 0;
  double density_final = 0;
  double density_inc = 0;
  double temperature = 0;
  double power = 0;
  ;
  double a_cst = 0;
  std::string potential_type = "";
  /* Whether or not to perform a reverse compression run */
  bool reverse_comp = false;
};

struct test_type {
  bool is_testing = false;
};

class md_options_interface {
 public:
  static int mdmain(std::string xml_file);

  static int get_switch(options_type& options, simulation_type& sim);

  static int load_constructor_options(constructor_type& constructor_options);

  static int load_additional_options(options_type& options);

  static int load_simulation_options(simulation_type& sim_options);

  static int load_test_options(test_type& test);
};
