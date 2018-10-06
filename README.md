# Molecular Dynamics Simulation

A program simulating a Molecular Dynamics (MD) fluid under the influence of an only repulsive, modified Lennard-Jones potential

![first eq](http://latex.codecogs.com/gif.latex?%5Cphi_%7Bij%7D%20%28r%29%20%3D%20%5Cvarepsilon%5Cbigg%28%5Cdfrac%7B%5Csigma%7D%7B%28r%5E%7B2%7D%20&plus;%20A%29%7D%5Cbigg%29%5E%5Cfrac%7Bn%7D%7B2%7D)


## Description
A program written to investigate the transition of a fluid from MD to SPH. The C++ files **MD.cpp** and **MD.h** contain the core of the fluid simulation while **FileLoadingPlotting.py** contains the analysis methods. The files **Source.cpp** and **Run.py** are used to execute the fluid simulation and the fluid analysis, repsectively.

## Getting Started
Firstly, the **FileIO** header file is required, which can be found [here](https://github.com/GiannisNikiteas/FileIO.git), clone it to the tools directory of this repository, like so:
```
# In the MD-simulation/tools directory
git clone https://github.com/GiannisNikiteas/FileIO.git
```

Then build the library and the test cases with, (change the compiler in **Makefile.variables** if needed):
```
make -j
```

### Path
* Change *path* in the constructor of **MD** to choose where the files will be saved.


### Run
Code to run the program can be found in **src/Source.cpp** and in the **examples** directory.
A minimal working example of how to use the **libmd.a** is shown below:

```
#include "MD.h"

int main()
{
  int steps = 5000; 	// number of iterations
  double a = 0.75; 	  // softness parameter
  int n = 6;  	      // potential strength
  double rho = 0.5; 	// fluid density
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
