# Molecular Dynamics Simulation

A program simulating a Molecular Dynamics (MD) fluid under the influence of an only repulsive, modified Lennard-Jones potential

![first eq](http://latex.codecogs.com/gif.latex?%5Cphi_%7Bij%7D%20%28r%29%20%3D%20%5Cvarepsilon%5Cbigg%28%5Cdfrac%7B%5Csigma%7D%7B%28r%5E%7B2%7D%20&plus;%20A%29%7D%5Cbigg%29%5E%5Cfrac%7Bn%7D%7B2%7D)



## Description
A program written to investigate the transition of a fluid from MD to SPH. The C++ files **MD.cpp** and **MD.h** contain the core of the fluid simulation while **FileLoadingPlotting.py** contains the analysis methods. The files **Source.cpp** and **Run.py** are used to execute the fluid simulation and the fluid analysis, repsectively.

## Getting Started
Download, **MD.cpp, MD.h, FileLoadingPlotting.py** and **Run.py**. **Source.cpp**, provides insight on how the simulation is supposed to be run.

### Path
* Change *path* in the constructor of **MD** to choose where the files will be saved.

* Chose that same *path* in the **Run.py**

### Run
Code to run the program can be found in **Source.cpp**. The example above, demonstrates the basic idea.

```
#include "MD.h"

int main()
{
  int steps = 5000;
  float A = 0.75;
  int n = 6;

  MD run(n, A, steps);
  run.Simulation();

  return 0;
}
```



## Components
The simulation is comprosed by the following:
* master: Fluid in a Non-Isothermal container

* Isothermal: Fluid in an Isothermal container

* Force_Interaction: A constant force in the x-direction is applied to the particle that is closest to the centre of the box

* Particle_Tracking: Writes the x,y,z position and velovities of all the particles for all the steps. Used in a combination with Python file **particles plot.py**

* Force_Interaction_Particle_Tracking: a combination of the ***Force_Interaction*** and ***Particle_Tracking*** branches

## Results
Snapshot of the fluid in 2D
![2Dsnap](https://github.com/GiannisNikiteas/MD-simulation/blob/master/snap_n6_a075.png?raw=true)

3D snapshot of the fluid
![3Dsnap](https://github.com/GiannisNikiteas/MD-simulation/blob/master/3d075.png?raw=true)

![2d-many-particles](https://github.com/GiannisNikiteas/MD-simulation/blob/master/snap_2744_n6_a4.png?raw=true)

<!-- 2D Animation of the fluid with time
![vid]() -->
