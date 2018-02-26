#! /bin/bash


icpc -std=c++17 -Ofast -use-intel-optimized-headers -m64 -mkl=sequential -fp-model=precise -ipo -shared-intel -ipp-link=static -qopt-subscript-in-range -qopenmp -fstack-protector -mtune=core-avx-i MD.h MD.cpp Source.cpp -o a.out
# -mgpu-asm-dump \
# -daal \
# -tbb \
