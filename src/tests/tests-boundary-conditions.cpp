#pragma once
#include <catch.hpp>

#define protected public
#include "../MD.h"

SCENARIO("Boundary conditions", "[bcs]") {
  GIVEN("Periodic boundaries") {
    MD obj;

    obj.options.Lx = obj.options.Ly = obj.options.Lz = 1.0;
    obj.options.N = 1;

    WHEN("x, y, z are smaller than 0") {
      vector_3d<double> r({-0.1}, {-0.2}, {-0.3});

      REQUIRE(r.size() == obj.options.N);

      INFO("vector before applying BCs: " << r);

      obj.periodic_boundary_conditions(r);

      INFO("vector after applying BCs: " << r);

      THEN("the positions are moved to the right side of the box") {
        REQUIRE(r.x[0] == Approx(0.9));
        REQUIRE(r.y[0] == Approx(0.8));
        REQUIRE(r.z[0] == Approx(0.7));
      }
    }
    WHEN("x, y, z are greater than the box length") {
      vector_3d<double> r({1.1}, {1.2}, {1.3});

      REQUIRE(r.size() == obj.options.N);

      INFO("vector before applying BCs: " << r);

      obj.periodic_boundary_conditions(r);

      INFO("vector after applying BCs: " << r);

      THEN("the positions are moved to the left side of the box") {
        REQUIRE(r.x[0] == Approx(0.1));
        REQUIRE(r.y[0] == Approx(0.2));
        REQUIRE(r.z[0] == Approx(0.3));
      }
    }
  }
  GIVEN("Reflective boundaries") {
    MD obj;

    obj.options.Lx = obj.options.Ly = obj.options.Lz = 1.0;
    obj.options.N = 1;

    WHEN("x, y, z are smaller than 0") {
      vector_3d<double> r({-0.1}, {-0.2}, {-0.3});
      vector_3d<double> v({-1}, {-1}, {-1});

      REQUIRE(r.size() == obj.options.N);
      REQUIRE(v.size() == r.size());

      INFO("vector before applying BCs: " << v);

      obj.reflective_boundary_conditions(r, v);

      INFO("vector after applying BCs: " << v);

      THEN("the positions are moved to the right side of the box") {
        REQUIRE(v.x[0] == Approx(1));
        REQUIRE(v.y[0] == Approx(1));
        REQUIRE(v.z[0] == Approx(1));
      }
    }
    WHEN("x, y, z are greater than the box length") {
      vector_3d<double> r({1.1}, {1.2}, {1.3});
      vector_3d<double> v({-1}, {-1}, {-1});

      REQUIRE(r.size() == obj.options.N);
      REQUIRE(v.size() == r.size());

      INFO("vector before applying BCs: " << v);

      obj.reflective_boundary_conditions(r, v);

      INFO("vector after applying BCs: " << v);

      THEN("the positions are moved to the left side of the box") {
        REQUIRE(v.x[0] == Approx(1));
        REQUIRE(v.y[0] == Approx(1));
        REQUIRE(v.z[0] == Approx(1));
      }
    }
  }
}