cmake_minimum_required(VERSION 2.6)

# Programming languages and system libraries
enable_language(Fortran)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Libraries
set(module_bvp       "bvp/bvp_la.f" "bvp/bvp_m.f90")
set(module_conductor "modules/module_precision.f90" "modules/module_spin.f90" "modules/module_state.f90" "modules/module_conductor.f90")

# Executables
add_executable(test "test.f90" ${module_conductor} ${module_bvp})
target_link_libraries(test ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
