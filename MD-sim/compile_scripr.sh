#! /bin/bash

icpc -std=c++17 -O3 -fast -pthread -mkl cluster -daal -parallel -fp-model consistent -prec-sqrt -mgpu-asm-dump -m64 -vec-guard-write MD.cpp stat_analysis.cpp Source.cpp -o a.out

./a.out
