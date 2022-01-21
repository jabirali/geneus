# This file is used by `cmake` to generate a `Makefile` for the project.
# It is recommended that you use the `configure` script to invoke `cmake`.

# Project metadata.
project(GENEUS Fortran)
cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

# Only one build type, compiler, and platform is currently supported.
set(CMAKE_BUILD_TYPE Debug)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(NOT Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message(WARNING "Only the GNU Fortran compiler is officially supported.")
endif()
if(NOT UNIX)
  message(WARNING "Only Unix-like operating systems are officially supported.")
endif()

# Customize the compile and linking options.
set(EXECUTABLE_OUTPUT_PATH "../bin/")
set(CMAKE_Fortran_FLAGS    "-march=x86-64 -Og -frecursive -fno-strict-aliasing -fno-stack-arrays -g -pg -fbacktrace")
set(CMAKE_EXE_LINKER_FLAGS "-march=x86-64 -Og -frecursive -fno-strict-aliasing -fno-stack-arrays -static -static-libgfortran")

# Customize the compile flags per source type.
file(GLOB lib77 RELATIVE ${CMAKE_SOURCE_DIR} "src/external/*.f")
foreach(src_name ${lib77})
  set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -std=legacy -ffixed-form)
endforeach()

file(GLOB lib90 RELATIVE ${CMAKE_SOURCE_DIR} "src/external/*.f90")
foreach(src_name ${lib90})
  set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -std=legacy -ffree-form)
endforeach()

file(GLOB lib08 RELATIVE ${CMAKE_SOURCE_DIR} "src/system/*.f90" "src/foundation/*.f90" "src/materials/*.f90")
file(GLOB bin08 RELATIVE ${CMAKE_SOURCE_DIR} "src/programs/*.f90")
foreach(src_name ${lib08} ${bin08})
  set_property(SOURCE ${src_name} APPEND PROPERTY COMPILE_OPTIONS -std=f2008 -ffree-form -fimplicit-none -fcheck=all,no-array-temps -pedantic)
endforeach()

# Define the compile targets: one library, many programs.
add_library(geneus STATIC ${lib77} ${lib90} ${lib08})
foreach(src_name ${bin08})
  get_filename_component( exe_name ${src_name} NAME_WE )
  add_executable( ${exe_name} ${src_name} )
  target_link_libraries( ${exe_name} geneus )
endforeach()

# Print the results.
message(STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "Fortran flags: " ${CMAKE_Fortran_FLAGS})
message(STATUS "Linker flags: " ${CMAKE_EXE_LINKER_FLAGS})
