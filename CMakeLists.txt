cmake_minimum_required(VERSION 3.0)

project(oasim-lib LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()


set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
# Flags for using and checking do concurrent instrinsic 
set(VECTOR_FLAGS 
              -ftree-vectorize 
              -fopt-info-vec 
              -std=f2008 
              -cpp
              )
set(CMAKE_VERBOSE_MAKEFILE ON)
# Set compilation flags depending on the compiler used, currently supports: GNU and Intel compilers.
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message(STATUS "Using GNU Fortran")
  # Added -pg for profiling in the debug mode only works with gnu compiler 
  # NOTE: will add vectorization flags for better peformance 
  add_compile_options(-Wall -Wextra -pedantic  -pg ${VECTOR_FLAGS})
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -fexpensive-optimizations -fomit-frame-pointer -ffast-math")
  # Added -pg for profiling in debug mode  
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message(STATUS "Using Intel Fortran")
  add_compile_options(-warn all) 
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -fomit-frame-pointer") # check compatibility with MPI for -fast
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -cpp -CB -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces -fpe0 -extend-source")
endif()

add_subdirectory(OASIMlib)
add_subdirectory(test)
add_subdirectory(benchmark) 

