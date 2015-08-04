# This file is used to generate a makefile for the project.
#
# Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Created: 2015-07-12
# Updated: 2015-07-28



###############################################################################
#                            PROJECT DEFINITION
###############################################################################

# General metadata
project(s4tran Fortran)
cmake_minimum_required(VERSION 2.6)

# Programming languages
enable_language(Fortran)

# Use static libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBRARIES OFF)

# Source files
file(GLOB src_ext77 "src/external/bvp_la.f")
file(GLOB src_ext90 "src/external/bvp_m.f" "src/external/pchip.f")
file(GLOB src_mod   "src/foundation/*.f" "src/materials/*.f" "src/algorithms/*.f")
file(GLOB src_bin   "src/tests/*.f" "src/programs/*.f")

# Executables

# Module tests
add_executable(test_foundation "src/tests/test_foundation.f" ${src_mod} ${src_ext77} ${src_ext90})
add_executable(test_materials "src/tests/test_materials.f" ${src_mod} ${src_ext77} ${src_ext90})
add_executable(test_critical "src/tests/test_critical.f" ${src_mod} ${src_ext77} ${src_ext90})


###############################################################################
#                            COMPILATION SETTINGS
###############################################################################

# Set the target directory
set(EXECUTABLE_OUTPUT_PATH "../bin")

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
  set(CMAKE_Fortran_FLAGS_RELEASE "-march=native -Ofast -flto -fno-stack-arrays -s")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-march=native -Og -g -pg -fbacktrace")
  set(CMAKE_EXE_LINKER_FLAGS      "-march=native -flto -fwhole-program -static")

  # Flags that apply to external files (imported from other projects)
  set_source_files_properties(${src_ext77} PROPERTIES COMPILE_FLAGS "-std=legacy -ffixed-form")
  set_source_files_properties(${src_ext90} PROPERTIES COMPILE_FLAGS "-std=legacy -ffree-form")

  # Flags that apply to internal files (produced as part of this project)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_source_files_properties(${src_mod} ${src_bin} PROPERTIES COMPILE_FLAGS "-std=f2008 -ffree-form -fimplicit-none -fcheck=all -Wconversion -pedantic")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(${src_mod} ${src_bin} PROPERTIES COMPILE_FLAGS "-std=f2008 -ffree-form")
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
