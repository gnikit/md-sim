#include "phase_transition.h"

void phase_transition::crystallisation(std::string SIMULATION_NAME,
                                       double DENSITY, double FINAL_DENSITY,
                                       double DENSITY_INC, double TEMPERATURE,
                                       double POWER, double A_CST,
                                       std::string pp_type) {
  /*
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
  // todo: many features do not work at this moment like particle tracking
  // todo: the POS << time_stamp stream will be a mess, same with RDF
  double current_rho = DENSITY;
  double old_box_length = 0;

  try {
    if (FINAL_DENSITY > DENSITY && DENSITY_INC > 0) {
      std::runtime_error(
          "You are attempting to compress the fluid:"
          "Final density has to be smaller than initial density");
    }
    if (FINAL_DENSITY < DENSITY && DENSITY_INC < 0) {
      std::runtime_error(
          "You are attempting to melt the crystal:"
          "Final density has to be larger than initial density");
    }
    if (DENSITY_INC > FINAL_DENSITY) {
      std::runtime_error(
          "Density increment has to be smaller than final density");
    }
  } catch (const std::exception &msg) {
    std::cerr << "Error: " << msg.what() << std::endl;
    exit(1);
  }

  // Number of compressions to occur
  size_t total_comp_steps =
      abs(ceil((FINAL_DENSITY - DENSITY) / DENSITY_INC)) + 1;

  for (size_t comp_step = 0; comp_step < total_comp_steps; comp_step++) {
    std::cout << "Runing MD::simulation " << comp_step + 1 << "/"
              << total_comp_steps << std::endl;

    simulation(SIMULATION_NAME, current_rho, TEMPERATURE, POWER, A_CST,
               pp_type);

    // Holds the box length of the previous simulation just run
    old_box_length = __L;

    // Density incrementation
    current_rho += DENSITY_INC;

    // Simulation updates old_box_length
    // the updated current_rho can generate the new box length
    // This value gets recalculated in the next Simulation
    __L = pow((__N / current_rho), 1.0 / 3.0);
    double box_length_ratio = __L / old_box_length;

    // Rescalling the positional vectors
    for (size_t i = 0; i < __N; ++i) {
      rx[i] *= box_length_ratio;
      ry[i] *= box_length_ratio;
      rz[i] *= box_length_ratio;
    }
    ++c_counter;
  }

  reset_values(true);
}

void phase_transition::run_backwards(std::string SIMULATION_NAME,
                                     double DENSITY, double FINAL_DENSITY,
                                     double DENSITY_INC, double TEMPERATURE,
                                     double POWER, double A_CST,
                                     std::string pp_type) {
  /*
   * //todo: implement in schema
   *
   */

  std::cout << "***************************\n"
               "** Starting Forward Run: **\n"
               "***************************\n"
            << std::endl;

  crystallisation("forward_" + SIMULATION_NAME, DENSITY, FINAL_DENSITY,
                  DENSITY_INC, TEMPERATURE, POWER, A_CST, pp_type);

  std::cout << "****************************\n"
               "** Finished Forward Run:  **\n"
               "** Starting Backward Run: **\n"
               "****************************\n"
            << std::endl;

  crystallisation("backward_" + SIMULATION_NAME, FINAL_DENSITY, DENSITY,
                  -DENSITY_INC, TEMPERATURE, POWER, A_CST, pp_type);

  std::cout << "*****************************\n"
               "** Finished Backward Run: **\n"
               "*****************************\n"
            << std::endl;
}