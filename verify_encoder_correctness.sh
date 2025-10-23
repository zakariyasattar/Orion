#!/bin/sh

# rm -rf CMakeCache.txt CMakeFiles
# cmake .
make -j8 linear_code_verifier
./linear_code_verifier