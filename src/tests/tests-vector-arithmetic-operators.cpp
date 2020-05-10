#pragma once
#include <catch.hpp>
#include "../vector_arithmetic_operators.h"

using namespace std;

SCENARIO("Testing std::vector<int> overloaded arithmetic operators",
         "[data-structures]") {
  GIVEN("testing operator+=") {
    WHEN("adding std::vector<int> to std::vector<int>") {
      vector<int> a{-1, 0, 1};
      vector<int> b{1, 0, -1};

      a += b;

      THEN("vectors have been added on a per element basis") {
        vector<int> ref{0, 0, 0};

        REQUIRE(a == ref);
      }
    }
    WHEN("adding <int> to std::vector<int>") {
      vector<int> a{0, 0, 0};
      int b = 1;

      a += b;

      THEN("a constant has been added to all elements in the vector") {
        vector<int> ref = {1, 1, 1};

        REQUIRE(a == ref);
      }
    }
  }
  GIVEN("testing operator+") {
    WHEN("adding std::vector<int> to std::vector<int>") {
      vector<int> a{-1, 0, 1};
      vector<int> b{1, 0, -1};

      vector<int> c = a + b;

      THEN("vectors have been added on a per element basis") {
        vector<int> ref{0, 0, 0};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
    WHEN("adding <int> to std::vector<int>") {
      vector<int> a{0, 0, 0};
      int b = 1;

      vector<int> c = a + b;

      THEN("a constant has been added to all elements in the vector") {
        vector<int> ref = {1, 1, 1};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
  }

  GIVEN("testing operator-=") {
    WHEN("subtracting std::vector<int> from std::vector<int>") {
      vector<int> a{1, 0, 1};
      vector<int> b{1, 0, 1};

      a -= b;

      THEN("vectors have been substracted on a per element basis") {
        vector<int> ref{0, 0, 0};

        REQUIRE(a == ref);
      }
    }
    WHEN("subtracting <int> to std::vector<int>") {
      vector<int> a{0, 0, 0};
      int b = -1;

      a -= b;

      THEN("a constant has been substracted from all elements in the vector") {
        vector<int> ref = {1, 1, 1};

        REQUIRE(a == ref);
      }
    }
  }
  GIVEN("testing operator-") {
    WHEN("subtracting std::vector<int> from std::vector<int>") {
      vector<int> a{1, 0, 1};
      vector<int> b{1, 0, 1};

      vector<int> c = a - b;

      THEN("vectors have been substracted on a per element basis") {
        vector<int> ref{0, 0, 0};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
    WHEN("subtracting <int> to std::vector<int>") {
      vector<int> a{0, 0, 0};
      int b = -1;

      vector<int> c = a - b;

      THEN("a constant has been substracted from all elements in the vector") {
        vector<int> ref = {1, 1, 1};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
  }

  GIVEN("testing operator*=") {
    WHEN("multiplying std::vector<int> to std::vector<int>") {
      vector<int> a{-1, 0, 1};
      vector<int> b{5, 0, -3};

      a *= b;

      THEN("vectors have been multiplied on a per element basis") {
        vector<int> ref{-5, 0, -3};

        REQUIRE(a == ref);
      }
    }
    WHEN("multiplying <int> to std::vector<int>") {
      vector<int> a{1, 2, 3};
      int b = 5;

      a *= b;

      THEN("a constant has been multiplied to all elements in the vector") {
        vector<int> ref = {5, 10, 15};

        REQUIRE(a == ref);
      }
    }
  }
  GIVEN("testing operator*") {
    WHEN("multiplying std::vector<int> to std::vector<int>") {
      vector<int> a{-1, 0, 1};
      vector<int> b{5, 0, -3};

      vector<int> c = a * b;

      THEN("vectors have been multiplied on a per element basis") {
        vector<int> ref{-5, 0, -3};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
    WHEN("multiplying <int> to std::vector<int>") {
      vector<int> a{1, 2, 3};
      int b = 5;

      vector<int> c = a * b;

      THEN("a constant has been multiplied to all elements in the vector") {
        vector<int> ref = {5, 10, 15};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
  }

  GIVEN("testing operator/=") {
    WHEN("dividing std::vector<int> with std::vector<int>") {
      vector<int> a{10, 11, 6};
      vector<int> b{2, 1, -2};

      a /= b;

      THEN("vectors have been divided on a per element basis") {
        vector<int> ref{5, 11, -3};

        REQUIRE(a == ref);
      }
    }
    WHEN("dividing std::vector<int> with <int>") {
      vector<int> a{15, 10, 5};
      int b = 5;

      a /= b;

      THEN("all elements in the vector have been divided by a constant") {
        vector<int> ref = {3, 2, 1};

        REQUIRE(a == ref);
      }
    }
  }
  GIVEN("testing operator/") {
    WHEN("dividing std::vector<int> with std::vector<int>") {
      vector<int> a{10, 11, 6};
      vector<int> b{2, 1, -2};

      vector<int> c = a / b;

      THEN("vectors have been divided on a per element basis") {
        vector<int> ref{5, 11, -3};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
    WHEN("dividing std::vector<int> with <int>") {
      vector<int> a{15, 10, 5};
      int b = 5;

      vector<int> c = a / b;

      THEN("all elements in the vector have been divided by a constant") {
        vector<int> ref = {3, 2, 1};

        REQUIRE(c == ref);
        REQUIRE(c != a);
      }
    }
  }
}