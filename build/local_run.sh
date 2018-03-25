#! /bin/bash


icpc -std=c++17 -O3 -Ofast -xHost -complex-limited-range -use-intel-optimized-headers -parallel -m64 -daal -tbb -mkl=parallel -fp-model=precise -ipo -shared-intel -ipp-link=static -qopt-subscript-in-range -qopenmp -fstack-protector -vec-guard-write -vecabi=cmdtarget -ax=SSE3 -mtune=core-avx2 -qopt-multi-version-aggressive -qopt-matmul -qopenmp-lib=compat -qopenmp-link=static -qopenmp-threadprivate=legacy -prec-sqrt -fast-transcendentals MD.h MD.cpp Source.cpp -o a.out

