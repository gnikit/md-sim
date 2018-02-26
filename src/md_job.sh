#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=4gb


# Load intel compiler
module load intel-suite/2017.1 
pwd
# Copy files to $TMPDIR
dir="$HOME"/MD/MD-simulation/src
cd "$dir"
cp vx.txt vy.txt vz.txt MD.h MD.cpp stat_analysis.h stat_analysis.cpp Source.cpp $TMPDIR

cd $TMPDIR
pwd
# Dir created in case of early termination
OUTDIR="$HOME"/MD/data
mkdir "$OUTDIR"

# Compile
icpc -std=c++17 -O3 -use-intel-optimized-headers -qopenmp -pthread -parallel MD.cpp stat_analysis.cpp Source.cpp -o run_me.out
# Run
./run_me.out
# Copy files from $TMPDIR to $WORK
cp * "$WORK"