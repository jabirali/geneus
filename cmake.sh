#!/bin/bash

cmake -DCMAKE_Fortran_COMPILER=/opt/gcc5/bin/gfortran -DCMAKE_BUILD_TYPE=Release -H. -Brelease
cmake -DCMAKE_Fortran_COMPILER=/opt/gcc5/bin/gfortran -DCMAKE_BUILD_TYPE=Debug   -H. -Bdebug
