git submodule update --init --recursive
rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake
cmake .
make linearPC_multi
./linearPC_multi 21
