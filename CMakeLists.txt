# This file is used to generate a makefile for the project.
#
# Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Created: 2015-07-12
# Updated: 2015-07-12

cmake_minimum_required(VERSION 2.6)

# Programming languages
enable_language(Fortran)

# System libraries
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

# Source files
file(GLOB src_bvp "bvp/*.f" "bvp/*.f90")
file(GLOB src_mod "modules/*.f90")

# Executables
add_executable(test "test.f90" ${src_mod} ${src_bvp})
target_link_libraries(test ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Global compiler flags
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message("-- Using Fortran compiler: " ${CMAKE_Fortran_COMPILER})
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fno-f2c")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -fno-f2c -g")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Using Fortran flags: " ${CMAKE_Fortran_FLAGS_DEBUG})
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    message("-- Using Fortran flags: " ${CMAKE_Fortran_FLAGS_RELEASE})
  endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
else(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message("-- Using Fortran compiler: " ${CMAKE_Fortran_COMPILER})
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Using Fortran flags: " ${CMAKE_Fortran_FLAGS_DEBUG})
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    message("-- Using Fortran flags: " ${CMAKE_Fortran_FLAGS_RELEASE})
  endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Local compiler flags
set_source_files_properties(${src_bvp} PROPERTIES COMPILE_FLAGS "-std=legacy")
set_source_files_properties(${src_mod} PROPERTIES COMPILE_FLAGS "-std=f2008")
