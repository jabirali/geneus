# This file is used to generate a makefile for the project.
#
# Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Created: 2015-07-12
# Updated: 2016-03-02



###############################################################################
#                            PROJECT DEFINITION
###############################################################################

# General metadata
project(s4tran Fortran)
cmake_minimum_required(VERSION 2.6)

# Programming languages
enable_language(Fortran)

# System libraries
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

# Source files
file(GLOB src_e77 RELATIVE ${CMAKE_SOURCE_DIR} "fortran/external/bvp_la.f")
file(GLOB src_e90 RELATIVE ${CMAKE_SOURCE_DIR} "fortran/external/bvp_m.f" "fortran/external/pchip.f" "fortran/external/fparser.f")
file(GLOB src_mod RELATIVE ${CMAKE_SOURCE_DIR} "fortran/system/*.f" "fortran/foundation/*.f" "fortran/materials/*.f")
file(GLOB src_exe RELATIVE ${CMAKE_SOURCE_DIR} "fortran/programs/*.f")

# Libraries
add_library(s4tran STATIC ${src_mod} ${src_e90} ${src_e77})

# Executables
foreach(src_name ${src_exe})
  get_filename_component( exe_name ${src_name} NAME_WE )
  add_executable( ${exe_name} ${src_name} )
  target_link_libraries( ${exe_name} s4tran ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
endforeach()



###############################################################################
#                            COMPILATION SETTINGS
###############################################################################

# Set the target directory
set(EXECUTABLE_OUTPUT_PATH "../bin/")

# Set the build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Using default build type: 'Release'.")
  set(CMAKE_BUILD_TYPE Release)
else()
  message(STATUS "Using custom build type: '" ${CMAKE_BUILD_TYPE} "'.")
endif()

# Check which compiler we use
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(Fortran_COMPILER_NAME MATCHES "gfortran.*" AND CMAKE_SYSTEM_NAME MATCHES "Windows")
  # Cross-compiling with GNU GFORTRAN for WINDOWS
  set(CMAKE_Fortran_FLAGS          "-flto -Ofast -fno-stack-arrays -s")
  set(CMAKE_EXE_LINKER_FLAGS       "-flto -Ofast -fwhole-program -static")
  set_source_files_properties(${src_e77} PROPERTIES COMPILE_FLAGS "-std=legacy -ffixed-form")
  set_source_files_properties(${src_e90} PROPERTIES COMPILE_FLAGS "-std=legacy -ffree-form")
  set_source_files_properties(${src_mod} PROPERTIES COMPILE_FLAGS "-std=f2008  -ffree-form -fimplicit-none -Wconversion")
  set_source_files_properties(${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008  -ffree-form -fimplicit-none -Wconversion")
elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*" AND CMAKE_BUILD_TYPE MATCHES "Release")
  # Compiling with GNU GFORTRAN in RELEASE mode
  set(CMAKE_AR                     "gcc-ar")
  SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
  SET(CMAKE_Fortran_ARCHIVE_FINISH  true)
  set(CMAKE_Fortran_FLAGS          "-march=native -flto -Ofast -fno-stack-arrays -s")
  set(CMAKE_EXE_LINKER_FLAGS       "-march=native -flto -Ofast -fwhole-program")
  set_source_files_properties(${src_e77} PROPERTIES COMPILE_FLAGS "-std=legacy -ffixed-form")
  set_source_files_properties(${src_e90} PROPERTIES COMPILE_FLAGS "-std=legacy -ffree-form")
  set_source_files_properties(${src_mod} PROPERTIES COMPILE_FLAGS "-std=f2008  -ffree-form -fimplicit-none -Wconversion")
  set_source_files_properties(${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008  -ffree-form -fimplicit-none -Wconversion")
elseif(Fortran_COMPILER_NAME MATCHES "gfortran.*" AND CMAKE_BUILD_TYPE MATCHES "Debug")
  # Compiling with GNU GFORTRAN in DEBUG mode
  set(CMAKE_AR                     "gcc-ar")
  SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
  SET(CMAKE_Fortran_ARCHIVE_FINISH  true)
  set(CMAKE_Fortran_FLAGS          "-march=native -Og -g -pg -fbacktrace")
  set(CMAKE_EXE_LINKER_FLAGS       "-march=native -Og -fwhole-program")
  set_source_files_properties(${src_e77} PROPERTIES COMPILE_FLAGS "-std=legacy -ffixed-form")
  set_source_files_properties(${src_e90} PROPERTIES COMPILE_FLAGS "-std=legacy -ffree-form")
  set_source_files_properties(${src_mod} PROPERTIES COMPILE_FLAGS "-std=f2008  -ffree-form -fimplicit-none -Wconversion -pedantic -fcheck=all")
  set_source_files_properties(${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008  -ffree-form -fimplicit-none -Wconversion -pedantic -fcheck=all")
elseif(Fortran_COMPILER_NAME MATCHES "ifort.*" AND CMAKE_BUILD_TYPE MATCHES "Release")
  # Compiling with INTEL IFORT in RELEASE mode
  set(CMAKE_Fortran_FLAGS    "-fast")
  set(CMAKE_EXE_LINKER_FLAGS "-fast")
  set_source_files_properties(${src_e77} PROPERTIES COMPILE_FLAGS "-stand none -fixed")
  set_source_files_properties(${src_e90} PROPERTIES COMPILE_FLAGS "-stand none -free")
  set_source_files_properties(${src_mod} PROPERTIES COMPILE_FLAGS "-stand f08  -free -warn declarations,interfaces,usage -diag-disable 5112,5268")
  set_source_files_properties(${src_exe} PROPERTIES COMPILE_FLAGS "-stand f08  -free -warn declarations,interfaces,usage -diag-disable 5112,5268")
elseif(Fortran_COMPILER_NAME MATCHES "ifort.*" AND CMAKE_BUILD_TYPE MATCHES "Debug")
  # Compiling with INTEL IFORT in DEBUG mode
  set(CMAKE_Fortran_FLAGS    "-xHost -g -traceback -ftrapuv")
  set(CMAKE_EXE_LINKER_FLAGS "-xHost")
  set_source_files_properties(${src_e77} PROPERTIES COMPILE_FLAGS "-stand none -fixed")
  set_source_files_properties(${src_e90} PROPERTIES COMPILE_FLAGS "-stand none -free")
  set_source_files_properties(${src_mod} PROPERTIES COMPILE_FLAGS "-stand f08  -free -warn declarations,interfaces,usage -check bounds,pointers,uninit,stack -diag-disable 5112,5268,10182")
  set_source_files_properties(${src_exe} PROPERTIES COMPILE_FLAGS "-stand f08  -free -warn declarations,interfaces,usage -check bounds,pointers,uninit,stack -diag-disable 5112,5268,10182")
else()
  # Compiling with an unknown compiler
  message(WARNING "The only compilers that have been tested are GFortran and IFort...")
endif()

# Print the results
message(STATUS "Using Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "Using Fortran flags: "    ${CMAKE_Fortran_FLAGS})
message(STATUS "Using linker flags: "     ${CMAKE_EXE_LINKER_FLAGS})
