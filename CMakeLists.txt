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

# Use static libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBRARIES OFF)

# System libraries
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

# Source files
file(GLOB src_ext "src/external/*.f" "src/external/*.f90")
file(GLOB src_mod "src/foundation/*.f90" "src/materials/*.f90" "src/algorithms/*.f90")
file(GLOB src_bin "src/tests/*.f90" "src/programs/*.f90")

# Executables

# Module tests
add_executable(test_foundation "src/tests/test_foundation.f90" ${src_mod} ${src_ext})
target_link_libraries(test_foundation ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(test_materials "src/tests/test_materials.f90" ${src_mod} ${src_ext})
target_link_libraries(test_materials ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(test_critical "src/tests/test_critical.f90" ${src_mod} ${src_ext})
target_link_libraries(test_critical ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


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
  set(CMAKE_Fortran_FLAGS_RELEASE "-march=native -Ofast -flto -fno-stack-arrays -s")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-march=native -Og -g -pg -fbacktrace")
  set(CMAKE_EXE_LINKER_FLAGS      "-march=native -flto -fwhole-program -static")

  # Flags that apply to external files (imported from other projects)
  set_source_files_properties(${src_ext} PROPERTIES COMPILE_FLAGS "-std=legacy")

  # Flags that apply to internal files (produced as part of this project)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_source_files_properties(${src_mod} ${src_bin} PROPERTIES COMPILE_FLAGS "-std=f2008 -Wconversion -pedantic -fimplicit-none -fcheck=all")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(${src_mod} ${src_bin} PROPERTIES COMPILE_FLAGS "-std=f2008")
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
