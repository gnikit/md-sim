#! /bin/bash

#PBS -l walltime=00:00:30
#PBS -l select=1:ncpus=4:mem=2gb


# Load intel compiler
module load intel-suite/2017.1 

# Loads initial velocities

cp "$HOME"/MD/MD-simulation/MD\ simulation/vx.txt $TEMP
cp "$HOME"/MD/MD-simulation/MD\ simulation/vy.txt $TEMP
cp "$HOME"/MD/MD-simulation/MD\ simulation/vy.txt $TEMP


icpc -std=c++17 -O3 MD.cpp stat_analysis.cpp Source.cpp -o run_me.out
./run_me.out
OUTDIR="$WORK"/data
mkdir "$OUTDIR"
# Copy files from temp to local file
cp -r "$TEMP"/data/ "$OUTDIR"