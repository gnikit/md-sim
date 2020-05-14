#pragma once
#include <catch.hpp>

#include "../mdmain.h"

using namespace Spud;
using namespace Catch::literals;

SCENARIO("MD XML interface", "[xml] [MD] [constructor]") {
  load_options("mdmain_test_options.xml");
  options_type opts;

  GIVEN("setup options") {
    md_options_interface::load_setup_options(opts);

    WHEN("loading mandatory setup options") {

      THEN("mandatory setup options have been loaded") {
        REQUIRE(opts.dimension == 3);
        REQUIRE(opts.steps == 100);
        REQUIRE_THAT(opts.particles, Catch::Equals<size_t>({2, 2, 2}));
        REQUIRE(opts.lattice == "FCC");
        REQUIRE(opts.iterative_method == "Verlet");
        REQUIRE(opts.bcs == "Periodic");
      }
    }

    WHEN("loading optional RDF setup options") {
      THEN("optional RDF setup have been loaded") {
        REQUIRE(opts.rdf_options.rdf_bins == 100);
        REQUIRE(opts.rdf_options.rdf_wait == 10);
      }
    }

    WHEN("loading optional box dimensions options") {
      THEN("optional box dimensions have been loaded") {
        REQUIRE(opts.Lx == 1.0_a);
        REQUIRE(opts.Ly == 1.5_a);
        REQUIRE(opts.Lz == 2.0_a);
      }
    }
  }

  GIVEN("io options") {
    md_options_interface::load_io_options(opts.io_options);

    WHEN("loading io options (sim name and output dir)") {
      THEN("(sim name and output dir) io options have been loaded") {
        REQUIRE(opts.io_options.simulation_name == "test_");
        REQUIRE(opts.io_options.dir == ".");
      }
    }

    WHEN("load true/false io options") {
      THEN("true/false io options have been loaded") {
        REQUIRE(opts.io_options.visualise == true);
        REQUIRE(opts.io_options.energies == true);
        REQUIRE(opts.io_options.pressure == true);
        REQUIRE(opts.io_options.msd == true);
        REQUIRE(opts.io_options.vaf == true);
        REQUIRE(opts.io_options.sf == true);
        REQUIRE(opts.io_options.rdf == true);
        REQUIRE(opts.io_options.position == true);
        REQUIRE(opts.io_options.compression_stats == true);
      }

      THEN("all true/false io options have been set to true") {
        REQUIRE_FALSE(opts.io_options.visualise == false);
        REQUIRE_FALSE(opts.io_options.energies == false);
        REQUIRE_FALSE(opts.io_options.pressure == false);
        REQUIRE_FALSE(opts.io_options.msd == false);
        REQUIRE_FALSE(opts.io_options.vaf == false);
        REQUIRE_FALSE(opts.io_options.sf == false);
        REQUIRE_FALSE(opts.io_options.rdf == false);
        REQUIRE_FALSE(opts.io_options.position == false);
        REQUIRE_FALSE(opts.io_options.compression_stats == false);
      }
    }
  }

  GIVEN("simulation options") {
    md_options_interface::load_simulation_options(opts);

    WHEN("loading simulation options") {
      THEN("simulation options have been loaded") {
        REQUIRE(opts.simulation_type == "NormalRun");
        REQUIRE(opts.potential_type == "BoundedInversePower");
        REQUIRE(opts.density == 0.5_a);
        REQUIRE(opts.target_temperature == 1.0_a);
        REQUIRE(opts.power == 6_a);
        REQUIRE(opts.a_cst == 1.0_a);
        REQUIRE(opts.dt == 0.1_a);
        REQUIRE(opts.normalise_dt_w_temp == true);
        REQUIRE(opts.cut_off == 2.0_a);
      }
    }
  }

  GIVEN("test options") {
    md_options_interface::load_test_options(opts.test_options);

    WHEN("loading test optins") {
      THEN("test options have been loaded") {
        REQUIRE(opts.test_options.is_testing == true);
      }
    }
  }
}