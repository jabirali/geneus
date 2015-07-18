# This file is used to generate a makefile for the project.
#
# Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Created: 2015-07-12
# Updated: 2015-07-16



###############################################################################
#                            PROJECT DEFINITION
###############################################################################

# General metadata
project(s4tran)
cmake_minimum_required(VERSION 2.6)

# Programming languages
enable_language(Fortran)

# System libraries
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

# Source files
file(GLOB src_bvp "external/*.f" "external/*.f90")
file(GLOB src_mod "modules/*.f90")
file(GLOB src_exe "tests/*.f90")

# Executables

# Module tests
add_executable(test_core "tests/test_core.f90" ${src_mod} ${src_bvp})
target_link_libraries(test_core ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(test_materials "tests/test_materials.f90" ${src_mod} ${src_bvp})
target_link_libraries(test_materials ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})



###############################################################################
#                            COMPILATION SETTINGS
###############################################################################

# Set the target directory
#set(EXECUTABLE_OUTPUT_PATH "bin")

# Set the build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Using default build type: 'Release'.")
  set(CMAKE_BUILD_TYPE Release)
else()
  message(STATUS "Using custom build type: '" ${CMAKE_BUILD_TYPE} "'.")
endif()

# Set the compiler flags
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # Flags that apply to all files
  set(CMAKE_Fortran_FLAGS_RELEASE "-march=native -Ofast -flto -fno-stack-arrays")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-march=native -Og -g -pg -fbacktrace")
  set(CMAKE_EXE_LINKER_FLAGS      "-march=native -flto -fwhole-program")

  # Flags that apply to external files (imported from other projects)
  set_source_files_properties(${src_bvp} PROPERTIES COMPILE_FLAGS "-std=legacy")

  # Flags that apply to internal files (produced as part of this project)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_source_files_properties(${src_mod} ${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008 -Wconversion -pedantic -fimplicit-none -fcheck=all")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(${src_mod} ${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008")
  endif()
endif()

# Print the results
message(STATUS "Using Fortran compiler: " ${CMAKE_Fortran_COMPILER})
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Using Fortran flags: " ${CMAKE_Fortran_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Using Fortran flags: " ${CMAKE_Fortran_FLAGS_RELEASE})
endif()
message(STATUS "Using linker flags: " ${CMAKE_EXE_LINKER_FLAGS})
