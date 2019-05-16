# Molecular Dynamics Simulation

[![Build Status](https://travis-ci.com/GNikit/MD-simulation.svg?branch=master)](https://travis-ci.com/GNikit/MD-simulation)

A program simulating Molecular Dynamics (MD) fluids, with the option to use
custom pair potentials. The potentials currently in the project include:

* Bounded Inverse Power (BIP)
* Gaussian Core Model (GCM)
* Exponential Pair Potential (EXP)
* Lennard-Jones Potential (LJ)

## Description

A program written to investigate the transition of a fluid from MD to SPH. The C++ files **MD.cpp** and **MD.h** contain the core of the fluid simulation while **tools/MD-Simulation-Data-Analysis** contains the analysis methods written in Python. The files **bin/md** executable in combination with the script under **tools/generate_xml_file.py** can be used to simulate a fluid with any initial parameters, by simply changing the xml file passed to **md**.

Examples are listed under the **examples** directory

## Getting Started

### Clone

You can clone the project with all its submodules with:

```bash
git clone --recurse-submodules https://github.com/GNikit/MD-simulation.git
```

### Build

Build the library, the executables and the test cases with:

```bash
make -j
```

The current configuration uses the Intel* compilers for C++, change the **Makefile.variables** to suit your compiler and compiler flags.

### Run

Code to run the program with an xml script can be found in **src/Source.cpp** and should be executed by:

```bash
./bin/md schemas/input_schema.xml
// ^--executable    ^-- input variables
```

Code without the use of the xml parser can be found in the **examples** directory.

### Use libmd in your project

A minimal working example of how to use the **libmd.a** is shown below, make sure to link and include the static library and the header file **MD.h** to your program.

```C++
#include "MD.h"

int main() {
  int steps = 5000;   // number of iterations
  double a = 0.75;    // softness parameter
  int n = 6;          // potential strength
  double rho = 0.5;   // fluid density
  double t = 0.5;     // fluid temperature
  // Files will be saved to
  std::string dir = "/Desired/Path/"; // The directory has to exist
  
  MD run(dir, steps); // More elaborate constructors also exist
  run.Simulation(rho, t, n, a, "LJ"); // LJ: Lennard-Jones pair potential
}
```
