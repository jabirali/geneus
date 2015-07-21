#!/bin/bash

rm -r debug release
cmake -DCMAKE_Fortran_COMPILER=/opt/gcc5/bin/gfortran -DCMAKE_C_COMPILER=/opt/gcc5/bin/gcc -DCMAKE_CXX_COMPILER=/opt/gcc5/bin/g++ -DCMAKE_BUILD_TYPE=Release -H. -Brelease
cmake -DCMAKE_Fortran_COMPILER=/opt/gcc5/bin/gfortran -DCMAKE_C_COMPILER=/opt/gcc5/bin/gcc -DCMAKE_CXX_COMPILER=/opt/gcc5/bin/g++ -DCMAKE_BUILD_TYPE=Debug   -H. -Bdebug
