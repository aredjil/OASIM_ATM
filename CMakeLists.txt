cmake_minimum_required(VERSION 3.0)
project(oasim-lib LANGUAGES Fortran)

option(USE_MPI "Build with MPI support" ON)
# Finding MPI if needed 
if(USE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif() 

# ---------------------------------------------------------
# Build type and module output
# ---------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(CMAKE_VERBOSE_MAKEFILE ON)

# ---------------------------------------------------------
# Profiling flags
# ---------------------------------------------------------
set(PROFILING_FLAGS "-pg")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# ---------------------------------------------------------
# Compiler-specific flags
# ---------------------------------------------------------
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    message(STATUS "Using GNU Fortran")
    add_compile_options(-Wall -Wextra -pedantic ${PROFILING_FLAGS})
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -fexpensive-optimizations -fomit-frame-pointer -ffast-math")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g ${PROFILING_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    message(STATUS "Using Intel Fortran")
    add_compile_options(-warn all ${PROFILING_FLAGS})
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -fomit-frame-pointer")
    set(CMAKE_Fortran_FLAGS_DEBUG   "--cpp -CB -check all -traceback -gen-interfaces -warn interfaces -fpe0 -extend-source ${PROFILING_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVHPC" OR
        "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
    message(STATUS "Using NVIDIA NVFortran")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast -acc -gpu=ccall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -acc -Minfo=accel -gpu=cc80")
    add_compile_options(-module ${CMAKE_Fortran_MODULE_DIRECTORY})
else()
    message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# ---------------------------------------------------------
# Subdirectories
# ---------------------------------------------------------
add_subdirectory(OASIMlib)
add_subdirectory(test)
add_subdirectory(benchmark)
