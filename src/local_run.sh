#! /bin/bash


icpc -std=c++17 -O3 -Ofast -use-intel-optimized-headers -m64 -daal -tbb -mkl=parallel -fp-model=precise -ipo -shared-intel -ipp-link=static -qopt-subscript-in-range -qopenmp -fstack-protector -mtune=core-avx-i -mgpu-asm-dump MD.h MD.cpp Source.cpp -o a.out

