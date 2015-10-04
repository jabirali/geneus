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
file(GLOB src_e77 RELATIVE ${CMAKE_SOURCE_DIR} "fortran/external/bvp_la.f")
file(GLOB src_e90 RELATIVE ${CMAKE_SOURCE_DIR} "fortran/external/bvp_m.f" "fortran/external/pchip.f")
file(GLOB src_mod RELATIVE ${CMAKE_SOURCE_DIR} "fortran/system/*.f" "fortran/foundation/*.f" "fortran/materials/*.f")
file(GLOB src_exe RELATIVE ${CMAKE_SOURCE_DIR} "fortran/programs/*.f")

# Libraries
add_library(s4tran STATIC ${src_mod} ${src_e90} ${src_e77})

# Executables
foreach(src_name ${src_exe})
  get_filename_component( exe_name ${src_name} NAME_WE )
  add_executable( ${exe_name} ${src_name} )
  target_link_libraries( ${exe_name} s4tran )
endforeach()



###############################################################################
#                            COMPILATION SETTINGS
###############################################################################

# Set the target directory
set(EXECUTABLE_OUTPUT_PATH "../programs")

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
  set(CMAKE_Fortran_FLAGS_DEBUG   "-march=native -Og -g -pg -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-march=native -Ofast -flto -fno-stack-arrays -s")
  set(CMAKE_EXE_LINKER_FLAGS      "-march=native -flto -fwhole-program -static")

  # Flags that apply to external files (imported from other projects)
  set_source_files_properties(${src_e77} PROPERTIES COMPILE_FLAGS "-std=legacy -ffixed-form")
  set_source_files_properties(${src_e90} PROPERTIES COMPILE_FLAGS "-std=legacy -ffree-form")

  # Flags that apply to internal files (produced as part of this project)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_source_files_properties(${src_mod} ${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008 -ffree-form -fimplicit-none -fcheck=all -Wconversion -pedantic")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(${src_mod} ${src_exe} PROPERTIES COMPILE_FLAGS "-std=f2008 -ffree-form")
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
