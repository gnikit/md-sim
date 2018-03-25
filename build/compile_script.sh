#! /bin/bash

icpc -std=c++17 -O3 -fast -pthread -I/src[...] -prec-sqrt -mgpu-asm-dump -m64 - MD.cpp stat_analysis.cpp Source.cpp -o a.out

# ./a.out
