#include <stdlib.h>
#include <any>
#include <map>
#include <string>
#include <thread>
#include "MD.h"
#include "phase_transition.h"
#include "spud"

class md_options_interface {
 public:
  static int mdmain(std::string xml_file);

  static int load_setup_options(options_type& options);

  static int load_simulation_options(options_type& options);

  static int load_test_options(test_options_type& test);
};
