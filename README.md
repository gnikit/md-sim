# Molecular Dynamics Simulation

A program simulating a Molecular Dynamics (MD) fluid under the influence of an only repulsive, modified Lennard-Jones potential

![first eq](http://latex.codecogs.com/gif.latex?%5Cphi_%7Bij%7D%20%28r%29%20%3D%20%5Cvarepsilon%5Cbigg%28%5Cdfrac%7B%5Csigma%7D%7B%28r%5E%7B2%7D%20&plus;%20A%29%7D%5Cbigg%29%5E%5Cfrac%7Bn%7D%7B2%7D)


## Description
A program written to investigate the transition of a fluid from MD to SPH. The C++ files **MD.cpp** and **MD.h** contain the core of the fluid simulation while **FileLoadingPlotting.py** contains the analysis methods. The files **Source.cpp** and **Run.py** are used to execute the fluid simulation and the fluid analysis, repsectively.

## Getting Started
Build the library and the test cases with:
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

## Results
Snapshot of the fluid in 2D
![2Dsnap](https://github.com/GiannisNikiteas/MD-simulation/blob/master/snap_n6_a075.png?raw=true)

3D snapshot of the fluid
![3Dsnap](https://github.com/GiannisNikiteas/MD-simulation/blob/master/3d075.png?raw=true)

![2d-many-particles](https://github.com/GiannisNikiteas/MD-simulation/blob/master/snap_2744_n6_a4.png?raw=true)

<!-- 2D Animation of the fluid with time
![vid]() -->
