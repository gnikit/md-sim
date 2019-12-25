#include <iostream>
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

  phase_transition(options_type &input_options);

  /**
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

  void set_compression_flag(bool is_compressing);
};