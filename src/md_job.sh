#! /bin/bash

#PBS -l walltime=00:00:30
#PBS -l select=1:ncpus=4:mem=2gb


# Load intel compiler
module load intel-suite/2017.1 
pwd
# Loads initial velocities
dir="$HOME"/MD/MD-simulation/src
cp "$dir"/{vx.txt, vy.txt, vz.txt, MD.h, MD.cpp, stat_analysis.h, stat_analysis.cpp, Source.cpp} $TEMP

cd $TEMP
pwd 

icpc -std=c++17 -O3 -use-intel-optimized-headers -qopenmp -pthread -parallel MD.cpp stat_analysis.cpp Source.cpp -o run_me.out
./run_me.out
OUTDIR="$WORK"/data
mkdir "$OUTDIR"
# Copy files from temp to local file
cp -r "$TEMP"/data/ "$OUTDIR"