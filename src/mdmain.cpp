#include "mdmain.h"

using namespace Spud;

int md_options_interface::mdmain(std::string xml_file) {
  // load the xml file in memory
  load_options(xml_file);
  constructor_type constructor;
  options_type options;
  simulation_type sim_options;
  test_type test_options;

  // load all options
  load_constructor_options(constructor);
  load_additional_options(options);
  load_simulation_options(sim_options);
  load_test_options(test_options);

  // determine which class to call
  if (sim_options.simulation_type == "NormalRun") {
    MD run(constructor.steps, constructor.particles, constructor.lattice);

    run.load_options(options.dir, options.track_particles, options.rdf_bins,
                     options.rdf_wait);

    if (test_options.is_testing) run.enable_testing(true);

    run.simulation(sim_options.simulation_name, sim_options.density,
                   sim_options.temperature, sim_options.power,
                   sim_options.a_cst, sim_options.potential_type);

  } else if (sim_options.simulation_type == "CompressionRun") {
    phase_transition run(constructor.steps, constructor.particles,
                         constructor.lattice);

    run.load_options(options.dir, options.track_particles, options.rdf_bins,
                     options.rdf_wait);

    if (test_options.is_testing) run.enable_testing(true);

    run.run_backwards(sim_options.simulation_name, sim_options.density,
                      sim_options.density_final, sim_options.density_inc,
                      sim_options.temperature, sim_options.power,
                      sim_options.a_cst, sim_options.potential_type);

  } else if (sim_options.simulation_type == "ReverseCompressionRun") {
    phase_transition run(constructor.steps, constructor.particles,
                         constructor.lattice);

    run.load_options(options.dir, options.track_particles, options.rdf_bins,
                     options.rdf_wait);

    if (test_options.is_testing) run.enable_testing(true);

    run.run_backwards(sim_options.simulation_name, sim_options.density,
                      sim_options.density_final, sim_options.density_inc,
                      sim_options.temperature, sim_options.power,
                      sim_options.a_cst, sim_options.potential_type);
  } else
    std::cerr << "Unrecognised simulation_name provided in xml" << std::endl;

  return 0;
}

int md_options_interface::load_constructor_options(
    constructor_type& constructor_options) {
  std::string path = "/constructor";
  OptionError error;

  int temp;
  std::vector<int> temp_v;

  error = get_option(path + "/steps", temp);
  assert(error == SPUD_NO_ERROR);
  constructor_options.steps = static_cast<size_t>(temp);

  error = get_option(path + "/particles", temp_v);
  assert(error == SPUD_NO_ERROR);
  constructor_options.particles.assign(temp_v.begin(), temp_v.end());

  error = get_option(path + "/lattice/name", constructor_options.lattice);
  assert(error == SPUD_NO_ERROR);

  return 0;
}

int md_options_interface::load_additional_options(options_type& options) {
  std::string path = "/options";
  OptionError error;
  int temp;

  if (have_option(path + "/output_dir")) {
    error = get_option(path + "/output_dir", options.dir);
    assert(error == SPUD_NO_ERROR);
  }

  if (have_option(path + "/rdf_bins")) {
    error = get_option(path + "/rdf_bins", temp);
    assert(error == SPUD_NO_ERROR);
    options.rdf_bins = static_cast<size_t>(temp);
  }

  if (have_option(path + "/track_particles")) options.track_particles = true;

  if (have_option(path + "/rdf_equilibrate")) {
    error = get_option(path + "/rdf_equilibrate", temp);
    assert(error == SPUD_NO_ERROR);
    options.rdf_wait = static_cast<size_t>(temp);
  }
  return 0;
}

int md_options_interface::load_simulation_options(
    simulation_type& sim_options) {
  std::string path = "/simulation_type";
  OptionError error;

  error = get_option(path + "/name", sim_options.simulation_type);
  assert(error == SPUD_NO_ERROR);

  error = get_option(path + "/simulation_name", sim_options.simulation_name);
  assert(error == SPUD_NO_ERROR);

  std::string pp_path = path + "/pair_potential";
  error = get_option(pp_path + "/name", sim_options.potential_type);
  assert(error == SPUD_NO_ERROR);

  error = get_option(pp_path + "/density", sim_options.density);
  assert(error == SPUD_NO_ERROR);

  error = get_option(pp_path + "/temperature", sim_options.temperature);
  assert(error == SPUD_NO_ERROR);

  if (have_option(pp_path + "/potential_strength")) {
    error =
        get_option(pp_path + "/potential_strength", sim_options.power);
    assert(error == SPUD_NO_ERROR);
  }
  if (have_option(pp_path + "/softening_parameter")) {
    error = get_option(pp_path + "/softening_parameter", sim_options.a_cst);
    assert(error == SPUD_NO_ERROR);
  }
  if (have_option(path + "/final_density")) {
    error = get_option(path + "/final_density", sim_options.density_final);
    assert(error == SPUD_NO_ERROR);
  }
  if (have_option(path + "/density_increment")) {
    error = get_option(path + "/density_increment", sim_options.density_inc);
    assert(error == SPUD_NO_ERROR);
  }
  return 0;
}

int md_options_interface::load_test_options(test_type& test) {
  std::string path = "/enable_testing";

  if (have_option(path)) test.is_testing = true;
  return 0;
}