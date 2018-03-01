#! /bin/bash

#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=4gb


# Load intel compiler
module load intel-suite/2017.1 

# Copy files to $TMPDIR
dir="$HOME"/MD/MD-simulation/src
cd "$dir"
cp vx.txt vy.txt vz.txt MD.h MD.cpp stat_analysis.h stat_analysis.cpp Source.cpp $TMPDIR

# IMPORTANT: SHELL DIR HAS TO BE SAME AS VX, VY, VZ
cd $TMPDIR
pwd
# Dir created in case of early termination
OUTDIR="$HOME"/MD/Data
mkdir "$OUTDIR"

# Compile
# icpc -std=c++17 -parallel -O3 -use-intel-optimized-headers MD.h MD.cpp Source.cpp -o a.out
icpc -fast -std=c++17 -parallel -use-intel-optimized-headers MD.cpp Source.cpp -o a.out
# Run
./a.out
# Copy files from $TMPDIR to $WORK

cp * "$OUTDIR"

# Display CPU model
lscpu
# Display intel compiler version
icpc --version
