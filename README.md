# Molecular Dynamics Simulation

A program simulating a Molecular Dynamics (MD) fluid under the influence of an only repulsive, modified Lennard-Jones potential

![first eq](http://latex.codecogs.com/gif.latex?%5Cphi_%7Bij%7D%20%28r%29%20%3D%20%5Cvarepsilon%5Cbigg%28%5Cdfrac%7B%5Csigma%7D%7B%28r%5E%7B2%7D%20&plus;%20A%29%7D%5Cbigg%29%5E%5Cfrac%7Bn%7D%7B2%7D)


## Description
A program written to investigate the transition of a fluid from MD to SPH. The C++ files **MD.cpp** and **MD.h** contain the core of the fluid simulation while **tools/MD-Simulation-Data-Analysis** contains the analysis methods written in Python. The files **bin/md** execuatable in combination with the script under **tools/generate_xml_file.py** can be used to simulate a fluid with any initial parameters, by simply changinh the xml file passed to **md**. 

Examples are listed under the **examples** directory 


## Getting Started

### Clone
You can clone the project with all its submodules with:
```
git clone --recurse-submodules -j4 https://github.com/GiannisNikiteas/MD-simulation.git
```

### Build
Build the library, the executables and the test cases with:
```
make -j
```
The current configuration uses the Intel* compilers for C++, change the **Makefile.variables** to suit your compiler and compiler flags.


### Run
Code to run the program with an xml script can be found in **src/Source.cpp** and should be executed by:
```
./bin/md schemas/input_schema.xml
// ^--executable    ^-- input variables
```
Code without the use of the xml parser can be found in the **examples** directory.

### Use libmd in your project
A minimal working example of how to use the **libmd.a** is shown below, make sure to link and include the static library and the header file **MD.h** to your program.

```
#include "MD.h"

int main() {
  int steps = 5000;   // number of iterations
  double a = 0.75;    // softness parameter
  int n = 6;          // potential strength
  double rho = 0.5;   // fluid density
  double t = 0.5;     // fluid temperature
  // Files will be saved to
  std::string dir = "/Desired/Path/"; // The directory has to exist
  
  MD run(dir, steps);
  run.Simulation(rho, t, n, a);
}
```
