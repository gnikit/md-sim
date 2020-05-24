#include <iostream>
#include <iterator>
#include "MD.h"

/* Load Intel math lib if available */
#if defined(__INTEL_COMPILER)
#include <mathimf.h> /* Intel Math library */
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

class phase_transition : public MD {
 public:
  using MD::MD;

 protected:
  std::ofstream compression_stats;

 public:
  phase_transition(options_type &input_options);

  /**
   * @brief
   * Compress the fluid to get the phase boundary for a specific temperature.
   *
   * ceil((FINAL_DENSITY - DENSITY) / DENSITY_INC) compressions of STEPS length
   * Performs repeated compresss of the fluid by periodically
   * incrementing the density of the fluid.
   * As a consequence the box length, the scaling factor and the
   * position vectors are also scaled in order to conserve the number
   * of particles in the box.
   *
   */
  void crystallisation(options_type &options);
  void crystallisation(std::string SIMULATION_NAME, double DENSITY,
                       double FINAL_DENSITY, double DENSITY_INC,
                       double TEMPERATURE, double POWER, double A_CST,
                       std::string pp_type);

  void two_way_compression(options_type &options);

  void detect_transition();

  /**
   * @brief Check the input options for any errors and either exit or
   * throw a warning, fix the problem and continue.
   * 
   * @param options 
   */
  void error_check_options(options_type &options);

  void set_compression_flag(bool is_compressing);
};