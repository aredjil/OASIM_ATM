cmake_minimum_required(VERSION 3.0)

project(oasim-lib LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

set(CMAKE_VERBOSE_MAKEFILE ON)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  add_compile_options(-Wall -Wextra -pedantic)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -fexpensive-optimizations -fomit-frame-pointer -ffast-math")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  add_compile_options(-warn all)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -fomit-frame-pointer") # check compatibility with MPI for -fast
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -cpp -CB -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces -fpe0 -extend-source")
endif()

add_subdirectory(OASIMlib)
add_subdirectory(test)
add_subdirectory(benchmark)

