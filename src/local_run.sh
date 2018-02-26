#! /bin/bash


icpc -std=c++17 -Ofast -pthread -use-intel-optimized-headers -m64 -mkl=parallel MD.h MD.cpp Source.cpp -o a.out
# -mgpu-asm-dump \
# -daal \
# -tbb \
