# This file is used to generate a makefile for the project.
#
# Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Created: 2015-07-12
# Updated: 2015-07-12

cmake_minimum_required(VERSION 2.6)

# Programming languages
enable_language(Fortran)
set (CMAKE_Fortran_FLAGS_RELEASE "-std=legacy -O3 -fno-f2c")
set (CMAKE_Fortran_FLAGS_DEBUG   "-std=legacy -O0 -fno-f2c -g")

# System libraries
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

# Define modules
set(module_bvp       "bvp/bvp_la.f" "bvp/bvp_m.f90")
set(module_core      "modules/module_precision.f90" "modules/module_spin.f90" "modules/module_state.f90")
set(module_conductor "modules/module_conductor.f90" ${module_core} ${module_bvp})

# Define executables
add_executable(test "test.f90" ${module_conductor})
target_link_libraries(test ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
