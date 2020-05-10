#pragma once
#include <catch.hpp>
#include <vector>

#include "../vector_3d.h"

SCENARIO("vector_3d can get populated", "[vector_3d], [data-structures]") {
  using Catch::Matchers::Approx;

  GIVEN("a double type") {
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    std::vector<double> z = {7, 8, 9};

    WHEN("calling the default constructor") {
      vector_3d<double> a;

      THEN("populate vectors with values") {
        a.x = x;
        a.y = y;
        a.z = z;

        // Check_that over required and specific Approx namespace to compare vec
        CHECK_THAT(a.x, Approx(x));
        CHECK_THAT(a.y, Approx(y));
        CHECK_THAT(a.z, Approx(z));
      }
    }
    WHEN("calling the initialising list constructor") {
      vector_3d<double> a(x, y, z);

      THEN("vectors have already populated within the constructors") {
        // Check_that over required and specific Approx namespace to compare vec
        CHECK_THAT(a.x, Approx(x));
        CHECK_THAT(a.y, Approx(y));
        CHECK_THAT(a.z, Approx(z));
      }
    }
  }
}

SCENARIO("vector_3d overloaded arithmetic operators",
         "[vector_3d], [data-structures]") {
  using Catch::Matchers::Approx;

  GIVEN("testing the operator+=") {

    WHEN("adding a vector_3d to vector_3d") {
      vector_3d<double> a({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
      vector_3d<double> b({1, 2, 3}, {4, 5, 6}, {7, 8, 9});

      a += b;

      THEN("vector_3d B added to vector_3d A") {
        std::vector<double> ref_x = {2, 4, 6};
        std::vector<double> ref_y = {8, 10, 12};
        std::vector<double> ref_z = {14, 16, 18};

        CHECK_THAT(a.x, Approx(ref_x));
        CHECK_THAT(a.y, Approx(ref_y));
        CHECK_THAT(a.z, Approx(ref_z));
      }
    }
    WHEN("adding a constant variable to a vector_3d") {
      vector_3d<double> a({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
      double b = 1;

      a += b;

      THEN("double B added to all elements of vector_3d A") {
        std::vector<double> ref_x = {2, 3, 4};
        std::vector<double> ref_y = {5, 6, 7};
        std::vector<double> ref_z = {8, 9, 10};

        CHECK_THAT(a.x, Approx(ref_x));
        CHECK_THAT(a.y, Approx(ref_y));
        CHECK_THAT(a.z, Approx(ref_z));
      }
    }
  }

  GIVEN("testing the operator-=") {

    WHEN("subtracting a vector_3d from another vector_3d") {
      vector_3d<double> a({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
      vector_3d<double> b({1, 2, 3}, {4, 5, 6}, {7, 8, 9});

      a -= b;

      THEN("vector_3d B subtracted from vector_3d A") {
        std::vector<double> ref_x = {0, 0, 0};
        std::vector<double> ref_y = {0, 0, 0};
        std::vector<double> ref_z = {0, 0, 0};

        CHECK_THAT(a.x, Approx(ref_x));
        CHECK_THAT(a.y, Approx(ref_y));
        CHECK_THAT(a.z, Approx(ref_z));
      }
    }
    WHEN("subtracting a double from a vector_3d") {
      vector_3d<double> a({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
      double b = 1;

      a -= b;

      THEN("double B subtracted from all elements of vector_3d A") {
        std::vector<double> ref_x = {0, 1, 2};
        std::vector<double> ref_y = {3, 4, 5};
        std::vector<double> ref_z = {6, 7, 8};

        CHECK_THAT(a.x, Approx(ref_x));
        CHECK_THAT(a.y, Approx(ref_y));
        CHECK_THAT(a.z, Approx(ref_z));
      }
    }
  }
}