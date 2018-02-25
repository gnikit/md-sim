#!usr/bin/bash
#PBS -l walltime=00:00:30
#PBS -l select=1:ncpus=4:mem=2gb


#Load intel compiler
module load intel-suite/2017.1 

cp $HOME/MD/MD-simulation/MD\ simulation/vx.txt $TEMP
cp $HOME/MD/MD-simulation/MD\ simulation/vy.txt $TEMP
cp $HOME/MD/MD-simulation/MD\ simulation/vy.txt $TEMP


