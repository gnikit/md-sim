#! /bin/bash


icpc -std=c++17 -O3 -parallel -pthread -use-intel-optimized-headers \
-prec-sqrt \
-qoffload=optional \
# -mgpu-asm-dump \
# -daal \
# -tbb \
-m64 \
MD.h MD.cpp stat_analysis.h stat_analysis.cpp Source.cpp \
-o a.out