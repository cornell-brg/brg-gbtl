#! /usr/bin/env bash

TOP=${PWD}

cd ${TOP}
rm -rf build-riscv
mkdir build-riscv
cd build-riscv
cmake ${TOP}/src -DCMAKE_CXX_COMPILER=clang++ -DARCH=RISCV -DPLATFORM=rvv_serial
make bfs_level_demo -i -j 8

cd ${TOP}
rm -rf build-x86
mkdir build-x86
cd build-x86
cmake ${TOP}/src -DCMAKE_CXX_COMPILER=clang++ -DPLATFORM=optimized_sequential
make bfs_level_ref_demo -i -j 8
