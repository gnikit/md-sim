#include "phase_transition.h"

phase_transition::phase_transition(options_type &input_options)
    : MD(input_options) {}

void phase_transition::crystallisation(options_type &options) {
  options.compression_options.compression = true;

  double current_rho = options.density;
  double old_box_length = 0;
  const std::string original_sim_name = options.io_options.simulation_name;

  error_check_options(options);

  /* Number of compressions to occur */
  size_t total_comp_steps =
      abs(ceil((options.compression_options.density_final - options.density) /
               options.compression_options.density_inc)) +
      1;

  /* Open a stream to output the compression stats summary */
  if (options.io_options.compression_stats) {
    /* The logger class is only used to generate the filename */
    stat_file temp_logger;
    std::string compression_stats_fname = temp_logger.file_naming(
        options.io_options.dir + "/" + options.io_options.simulation_name +
            "Compression_statistics",
        options.steps, options.N, options.density, options.target_temperature,
        options.power, options.a_cst);

    /* Open file stream */
    compression_stats.open(compression_stats_fname,
                           std::ios::out | std::ios::trunc);

    /* Create the header in the order that get_run_stats generates them
       @note This has been hard coded */
    std::string header = "# compress_step,rho";
    if (options.io_options.msd)
      header += ",min(MSD),max(MSD),av(MSD),l2norm(MSD),rms(MSD)";
    if (options.io_options.vaf)
      header += ",min(VAF),max(VAF),av(VAF),l2norm(VAF),rms(VAF)";
    if (options.io_options.sf) {
      header += ",min(SFx),max(SFx),av(SFx),l2norm(SFx),rms(SFx)";
      header += ",min(SFy),max(SFy),av(SFy),l2norm(SFy),rms(SFy)";
      header += ",min(SFz),max(SFz),av(SFz),l2norm(SFz),rms(SFz)";
    }
    if (options.io_options.energies)
      header += ",min(U),max(U),av(U),l2norm(U),rms(U)";
    if (options.io_options.pressure)
      header += ",min(Pc),max(Pc),av(Pc),l2norm(Pc),rms(Pc)";

    /* Write the header */
    compression_stats << header << std::endl;
  }

  for (size_t comp_step = 0; comp_step < total_comp_steps; comp_step++) {
    std::cout << "Runing MD::simulation " << comp_step + 1 << "/"
              << total_comp_steps << std::endl;

    std::string prefix = original_sim_name + std::to_string(comp_step) + "_";

    options.io_options.simulation_name = prefix;
    options.density = current_rho;
    MD::set_options(options);
    simulation();

    /* Write a summary of the stats of the simulation run */
    if (options.io_options.compression_stats) {
      std::vector<double> run_data = get_run_stats();
      stat_file temp_logger;
      temp_logger.write_file_line(run_stats, compression_stats, comp_step);
    }

    /* Holds the box length of the previous simulation just run */
    old_box_length = options.L;

    /* Density incrementation */
    current_rho += options.compression_options.density_inc;

    /* Simulation updates old_box_length
       the updated current_rho can generate the new box length
       This value gets recalculated in the next Simulation */
    options.L = pow((options.N / current_rho), 1.0 / 3.0);
    options.Lx = options.Ly = options.Lz = options.L;
    double box_length_ratio = options.L / old_box_length;
    options.cut_off *= box_length_ratio;

    /* Rescalling the positional vectors */
    for (size_t i = 0; i < options.N; ++i) {
      r.x[i] *= box_length_ratio;
      r.y[i] *= box_length_ratio;
      r.z[i] *= box_length_ratio;
    }
    ++options.compression_options.compress_count;
  }

  if (options.io_options.compression_stats) compression_stats.close();
  reset_values(true);
}

void phase_transition::crystallisation(std::string SIMULATION_NAME,
                                       double DENSITY, double FINAL_DENSITY,
                                       double DENSITY_INC, double TEMPERATURE,
                                       double POWER, double A_CST,
                                       std::string pp_type) {
  options.io_options.simulation_name = SIMULATION_NAME;
  options.density = DENSITY;
  options.compression_options.density_final = FINAL_DENSITY;
  options.compression_options.density_inc = DENSITY_INC;
  options.target_temperature = TEMPERATURE;
  options.power = POWER;
  options.a_cst = A_CST;
  options.potential_type = pp_type;

  crystallisation(options);
}

void phase_transition::two_way_compression(options_type &options) {
  /*
   * //todo: implement in schema
   *
   */

  std::cout << "***************************\n"
               "** Starting Forward Run: **\n"
               "***************************\n"
            << std::endl;

  crystallisation("forward_" + options.io_options.simulation_name,
                  options.density, options.compression_options.density_final,
                  options.compression_options.density_inc,
                  options.target_temperature, options.power, options.a_cst,
                  options.potential_type);

  std::cout << "****************************\n"
               "** Finished Forward Run:  **\n"
               "** Starting Backward Run: **\n"
               "****************************\n"
            << std::endl;

  crystallisation("backward_" + options.io_options.simulation_name,
                  options.compression_options.density_final, options.density,
                  -options.compression_options.density_inc,
                  options.target_temperature, options.power, options.a_cst,
                  options.potential_type);

  std::cout << "*****************************\n"
               "** Finished Backward Run: **\n"
               "*****************************\n"
            << std::endl;
}

void phase_transition::error_check_options(options_type &options) {
  /**************************** WARNINGS **************************************/
  try {
    if (options.compression_options.density_final > options.density) {
      std::cout << "Starting fluid compression" << std::endl;

      if (options.compression_options.density_inc < 0) {
        std::string warning =
            "Compressing the fluid requires a positive density increment to be "
            "passed. The input increment: " +
            std::to_string(options.compression_options.density_inc) +
            " will be replaced by: " +
            std::to_string(-options.compression_options.density_inc);

        // Change the sign of the density increment
        options.compression_options.density_inc =
            -options.compression_options.density_inc;

        throw std::runtime_error(warning);
      }
    } else if (options.compression_options.density_final < options.density) {
      std::cout << "Starting fluid decompression/ crystal melting" << std::endl;

      if (options.compression_options.density_inc > 0) {
        std::string warning =
            "Decompressing the fluid requires a negative density increment to "
            "be passed. The input increment: " +
            std::to_string(options.compression_options.density_inc) +
            " will be replaced by: " +
            std::to_string(-options.compression_options.density_inc);

        // Change the sign of the density increment
        options.compression_options.density_inc =
            -options.compression_options.density_inc;

        throw std::runtime_error(warning);
      }
    }
  } catch (const std::exception &msg) {
    std::cerr << "Warning: " << msg.what() << std::endl;
  }

  /****************************** ERRORS **************************************/
  try {
    if (options.compression_options.density_inc >
        options.compression_options.density_final) {
      throw std::runtime_error(
          "Density increment has to be smaller than final density");
    }
  } catch (const std::exception &msg) {
    std::cerr << "Error: " << msg.what() << std::endl;
    exit(-1);
  }

  /****************************************************************************/
}

void phase_transition::set_compression_flag(bool is_compressing) {
  options.compression_options.compression = is_compressing;
}
