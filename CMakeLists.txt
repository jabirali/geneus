# This file is used by `cmake` to generate a `Makefile` for the project.

# General metadata
project(GENEUS Fortran)
cmake_minimum_required(VERSION 3.0)

# Programming languages
enable_language(Fortran)

# Group source files by their Fortran version
file(GLOB lib77 RELATIVE ${CMAKE_SOURCE_DIR} "src/external/*.f")
file(GLOB lib90 RELATIVE ${CMAKE_SOURCE_DIR} "src/external/*.f90")
file(GLOB lib08 RELATIVE ${CMAKE_SOURCE_DIR} "src/system/*.f90" "src/foundation/*.f90" "src/materials/*.f90")
file(GLOB bin08 RELATIVE ${CMAKE_SOURCE_DIR} "src/programs/*.f90")

# Define the relevant compilation targets
add_library(geneus STATIC ${lib77} ${lib90} ${lib08})
foreach(src_name ${bin08})
  get_filename_component( exe_name ${src_name} NAME_WE )
  add_executable( ${exe_name} ${src_name} )
  target_link_libraries( ${exe_name} geneus )
endforeach()

# Set the target directory
set(EXECUTABLE_OUTPUT_PATH "../bin/")

# Only one build type, compiler, and platform is currently supported
set(CMAKE_BUILD_TYPE Debug)

if(NOT UNIX)
  message(WARNING "Only Unix-like operating systems are officially supported.")
endif()

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(NOT Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message(WARNING "Only the GNU Fortran compiler is officially supported.")
endif()

# Leniency for external code
foreach(src_name ${lib77})
  set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -std=legacy -ffixed-form)
endforeach()
foreach(src_name ${lib90})
  set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -std=legacy -ffree-form)
endforeach()

# Extra warnings for my code
foreach(src_name ${lib08} ${bin08})
  set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -std=f2008 -ffree-form -fimplicit-none -Wconversion)
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -pedantic -fcheck=all,no-array-temps)
  endif()
endforeach()

# Customize build flags
set(CMAKE_Fortran_FLAGS    "-march=native -Os -g -pg -fbacktrace")
set(CMAKE_EXE_LINKER_FLAGS "-march=native -Os -fwhole-program -static")

# Print the results
message(STATUS "Using Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "Using Fortran flags: "    ${CMAKE_Fortran_FLAGS})
message(STATUS "Using linker flags: "     ${CMAKE_EXE_LINKER_FLAGS})
