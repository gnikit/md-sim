#pragma once
#include <catch.hpp>

#define protected public
#include "../MD.h"

using Catch::Matchers::Approx;

SCENARIO("Lattice conditions", "[lattice]") {
  GIVEN("Simple Cubic lattice (SC)") {
    MD obj;

    obj.options.Lx = obj.options.Ly = obj.options.Lz = 1.0;
    obj.options.Nx = obj.options.Ny = obj.options.Nz = 2.0;
    obj.options.N = obj.options.Nx * obj.options.Ny * obj.options.Nz;

    vector_3d<double> r;
    r.reserve(obj.options.N);
    const std::string lattice = "SC";
    obj.choose_lattice_formation(lattice, r);

    THEN("positions have been arranged as a SC lattice") {
      INFO("positions are:\n" << r);

      vector_3d<double> ref({0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75},
                            {0.25, 0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75},
                            {0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75});

      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
    }
  }
  GIVEN("Face Centred Cubic lattice (FCC)") {
    MD obj;

    obj.options.Lx = obj.options.Ly = obj.options.Lz = 1.0;
    obj.options.Nx = obj.options.Ny = obj.options.Nz = 2.0;
    obj.options.N = obj.options.Nx * obj.options.Ny * obj.options.Nz * 4;

    vector_3d<double> r;
    r.reserve(obj.options.N);
    const std::string lattice = "FCC";
    obj.choose_lattice_formation(lattice, r);

    THEN("positions have been arranged as a FCC lattice") {
      INFO("positions are:\n" << r);

      vector_3d<double> ref(
          {0.125, 0.125, 0.125, 0.125, 0.625, 0.625, 0.625, 0.625,
           0.375, 0.375, 0.375, 0.375, 0.875, 0.875, 0.875, 0.875,
           0.375, 0.375, 0.375, 0.375, 0.875, 0.875, 0.875, 0.875,
           0.125, 0.125, 0.125, 0.125, 0.625, 0.625, 0.625, 0.625},
          {0.125, 0.125, 0.625, 0.625, 0.125, 0.125, 0.625, 0.625,
           0.375, 0.375, 0.875, 0.875, 0.375, 0.375, 0.875, 0.875,
           0.125, 0.125, 0.625, 0.625, 0.125, 0.125, 0.625, 0.625,
           0.375, 0.375, 0.875, 0.875, 0.375, 0.375, 0.875, 0.875},
          {0.125, 0.625, 0.125, 0.625, 0.125, 0.625, 0.125, 0.625,
           0.125, 0.625, 0.125, 0.625, 0.125, 0.625, 0.125, 0.625,
           0.375, 0.875, 0.375, 0.875, 0.375, 0.875, 0.375, 0.875,
           0.375, 0.875, 0.375, 0.875, 0.375, 0.875, 0.375, 0.875});

      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
    }
  }
  GIVEN("Body Centred Cubic lattice (BCC)") {
    MD obj;

    obj.options.Lx = obj.options.Ly = obj.options.Lz = 1.0;
    obj.options.Nx = obj.options.Ny = obj.options.Nz = 2.0;
    obj.options.N = obj.options.Nx * obj.options.Ny * obj.options.Nz * 2;

    vector_3d<double> r;
    r.reserve(obj.options.N);
    const std::string lattice = "BCC";
    obj.choose_lattice_formation(lattice, r);

    THEN("positions have been arranged as a BCC lattice") {
      INFO("positions are:\n" << r);

      vector_3d<double> ref(
          {0.125, 0.125, 0.125, 0.125, 0.625, 0.625, 0.625, 0.625, 0.375, 0.375,
           0.375, 0.375, 0.875, 0.875, 0.875, 0.875},
          {0.125, 0.125, 0.625, 0.625, 0.125, 0.125, 0.625, 0.625, 0.375, 0.375,
           0.875, 0.875, 0.375, 0.375, 0.875, 0.875},
          {0.125, 0.625, 0.125, 0.625, 0.125, 0.625, 0.125, 0.625, 0.375, 0.875,
           0.375, 0.875, 0.375, 0.875, 0.375, 0.875});

      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
    }
  }
  GIVEN("Random lattice option") {
    MD obj;

    obj.options.Lx = obj.options.Ly = obj.options.Lz = 1.0;
    obj.options.Nx = obj.options.Ny = obj.options.Nz = 2.0;
    obj.options.N = obj.options.Nx * obj.options.Ny * obj.options.Nz;
    obj.options.test_options.is_testing = true;

    vector_3d<double> r;
    r.reserve(obj.options.N);
    const std::string lattice = "RANDOM";
    obj.choose_lattice_formation(lattice, r);

    THEN("positions have been arranged as a RANDOM (uniform dist) manner") {
      INFO("positions are:\n" << r);

      vector_3d<double> ref({0.102618, 0.750331, 0.259876, 0.233401, 0.397476,
                             0.73461, 0.767412, 0.847782},
                            {0.211333, 0.335968, 0.281168, 0.767064, 0.316441,
                             0.731596, 0.953231, 0.349762},
                            {0.128977, 0.244343, 0.462842, 0.819957, 0.155512,
                             0.857859, 0.290974, 0.923897});

      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
      CHECK_THAT(r.x, Approx(ref.x));
    }
  }
}