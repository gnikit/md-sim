#! /bin/bash

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=4gb


# Load intel compiler
module load intel-suite/2017.1 

# Copy files to $TMPDIR
dir="$HOME"/MD/MD-simulation/src
cd "$dir"
cp vx.txt vy.txt vz.txt MD.h MD.cpp stat_analysis.h stat_analysis.cpp Source.cpp $TMPDIR

# IMPORTANT: SHELL DIR HAS TO BE SAME AS VX, VY, VZ FILES
cd $TMPDIR
pwd
# Dir created in case of early termination
OUTDIR="$HOME"/MD/data
mkdir "$OUTDIR"

# Compile
icpc -std=c++17 -O3 -Ofast -xHost -complex-limited-range -use-intel-optimized-headers -parallel -m64 -daal -tbb -mkl=parallel -fp-model=precise -ipo -shared-intel -ipp-link=static -qopt-subscript-in-range -qopenmp -fstack-protector -vec-guard-write -vecabi=cmdtarget -ax=SSE3 -mtune=core-avx2 -qopt-multi-version-aggressive -qopt-matmul -qopenmp-lib=compat -qopenmp-link=static -qopenmp-threadprivate=legacy -prec-sqrt -fast-transcendentals MD.h MD.cpp Source.cpp -o a.out
# Run
./a.out
# Copy files from $TMPDIR to $WORK

cp * "$OUTDIR"