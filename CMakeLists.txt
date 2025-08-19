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

  # Detect compiler
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  message(STATUS "Using GNU Fortran")
  add_compile_options(-Wall -Wextra -pedantic -pg ${VECTOR_FLAGS})
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -fexpensive-optimizations -fomit-frame-pointer -ffast-math")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message(STATUS "Using Intel Fortran")
  add_compile_options(-warn all)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -fomit-frame-pointer")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -cpp -CB -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces -fpe0 -extend-source")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVHPC" OR
        "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
  message(STATUS "Using NVIDIA NVFortran")
  # Enable OpenACC (you can adjust flags depending on your GPU target)
  set(CMAKE_Fortran_FLAGS_RELEASE "-fast -acc -gpu=ccall")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g -acc -Minfo=accel")
  # Optionally add module output directory
  add_compile_options(-module ${CMAKE_Fortran_MODULE_DIRECTORY})
else()
  message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

  add_subdirectory(OASIMlib)
  add_subdirectory(test)
  add_subdirectory(benchmark) 

