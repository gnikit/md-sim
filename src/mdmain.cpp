#include "mdmain.h"

using namespace Spud;

int md_options_interface::mdmain(std::string xml_file) {
  /* load the xml file in memory */
  load_options(xml_file);
  options_type opts;

  /* load all opts */
  load_setup_options(opts);
  load_io_options(opts.io_options);
  load_simulation_options(opts);
  load_test_options(opts.test_options);

  /* determine which class to call */
  if (opts.simulation_type == "NormalRun") {
    MD run(opts);

    run.simulation();

  } else if (opts.simulation_type == "CompressionRun") {
    phase_transition run(opts);

    run.crystallisation(opts);

  } else if (opts.simulation_type == "ReverseCompressionRun") {
    phase_transition run(opts);

    if (opts.test_options.is_testing) run.enable_testing(true);

    run.two_way_compression(opts);
  } else
    std::cerr << "Unrecognised simulation_name provided in xml" << std::endl;

  return 0;
}

int md_options_interface::load_setup_options(options_type& options) {
  std::string path = "/setup";
  OptionError error;

  int temp;
  std::vector<int> temp_v;

  error = get_option(path + "/steps", temp);
  assert(error == SPUD_NO_ERROR);
  options.steps = static_cast<size_t>(temp);

  error = get_option(path + "/particles", temp_v);
  assert(error == SPUD_NO_ERROR);
  options.particles.assign(temp_v.begin(), temp_v.end());

  error = get_option(path + "/lattice/name", options.lattice);
  assert(error == SPUD_NO_ERROR);

  error = get_option(path + "/iterative_method/name", options.iterative_method);
  assert(error == SPUD_NO_ERROR);

  error = get_option(path + "/boundary_conditions/name", options.bcs);
  assert(error == SPUD_NO_ERROR);

  if (have_option(path + "/rdf_bins")) {
    error = get_option(path + "/rdf_bins", temp);
    assert(error == SPUD_NO_ERROR);
    options.rdf_options.rdf_bins = static_cast<size_t>(temp);
  }

  if (have_option(path + "/rdf_equilibrate")) {
    error = get_option(path + "/rdf_equilibrate", temp);
    assert(error == SPUD_NO_ERROR);
    options.rdf_options.rdf_wait = static_cast<size_t>(temp);
  }

  if (have_option(path + "/track_particles"))
    options.io_options.visualise = true;

  return 0;
}

int md_options_interface::load_io_options(io_options_type& io) {
  std::string path = "/io";
  std::string temp;
  OptionError error;

  if (have_option(path + "/simulation_name")) {
    error = get_option(path + "/simulation_name", io.simulation_name);
    assert(error == SPUD_NO_ERROR);
  }

  if (have_option(path + "/output_dir")) {
    error = get_option(path + "/output_dir", io.dir);
    assert(error == SPUD_NO_ERROR);
  }

  error = get_option(path + "/track_particles/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.visualise = true;
  else
    io.visualise = false;

  error = get_option(path + "/energies/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.energies = true;
  else
    io.energies = false;

  error = get_option(path + "/pressure/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.pressure = true;
  else
    io.pressure = false;

  error = get_option(path + "/msd/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.msd = true;
  else
    io.msd = false;

  error = get_option(path + "/vaf/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.vaf = true;
  else
    io.vaf = false;

  error = get_option(path + "/sf/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.sf = true;
  else
    io.sf = false;

  error = get_option(path + "/rdf/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.rdf = true;
  else
    io.rdf = false;

  error = get_option(path + "/positions/name", temp);
  assert(error == SPUD_NO_ERROR);
  if (temp == "true")
    io.position = true;
  else
    io.position = false;

  if (have_option(path + "/compression_summary_stats")) {
    error = get_option(path + "/compression_summary_stats/name", temp);
    assert(error == SPUD_NO_ERROR);
    if (temp == "true")
      io.compression_stats = true;
    else
      io.compression_stats = false;
  }

  return 0;
}

int md_options_interface::load_simulation_options(options_type& opts) {
  std::string path = "/simulation_type";
  OptionError error;

  error = get_option(path + "/name", opts.simulation_type);
  assert(error == SPUD_NO_ERROR);

  std::string pp_path = path + "/pair_potential";
  error = get_option(pp_path + "/name", opts.potential_type);
  assert(error == SPUD_NO_ERROR);

  error = get_option(pp_path + "/density", opts.density);
  assert(error == SPUD_NO_ERROR);

  error = get_option(pp_path + "/temperature", opts.target_temperature);
  assert(error == SPUD_NO_ERROR);

  if (have_option(pp_path + "/potential_strength")) {
    error = get_option(pp_path + "/potential_strength", opts.power);
    assert(error == SPUD_NO_ERROR);
  }

  if (have_option(pp_path + "/softening_parameter")) {
    error = get_option(pp_path + "/softening_parameter", opts.a_cst);
    assert(error == SPUD_NO_ERROR);
  }

  if (have_option(pp_path + "/dt")) {
    error = get_option(pp_path + "/dt", opts.dt);
    assert(error == SPUD_NO_ERROR);

    if (have_option(pp_path + "/dt/normalise_with_temperature"))
      opts.normalise_dt_w_temp = true;
  } else
    opts.normalise_dt_w_temp = false;

  if (have_option(path + "/final_density")) {
    error = get_option(path + "/final_density",
                       opts.compression_options.density_final);
    assert(error == SPUD_NO_ERROR);
  }

  if (have_option(pp_path + "/cut_off")) {
    error = get_option(pp_path + "/cut_off", opts.cut_off);
    assert(error = SPUD_NO_ERROR);
  }

  if (have_option(path + "/density_increment")) {
    error = get_option(path + "/density_increment",
                       opts.compression_options.density_inc);
    assert(error == SPUD_NO_ERROR);
  }

  return 0;
}

int md_options_interface::load_test_options(test_options_type& test) {
  std::string path = "/enable_testing";

  if (have_option(path)) test.is_testing = true;

  return 0;
}