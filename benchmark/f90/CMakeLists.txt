cmake_minimum_required(VERSION 3.0)
option(ENABLE_PROFILING "Enable profiling with gprof" OFF)
if(ENABLE_PROFILING)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

execute_process(COMMAND nf-config --cflags
    OUTPUT_VARIABLE NETCDF_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get NetCDF linker flags
execute_process(COMMAND nf-config --flibs
    OUTPUT_VARIABLE NETCDF_FLIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
add_executable(benchmark benchmark.f90)

set_target_properties(benchmark PROPERTIES SUFFIX ".x")

separate_arguments(NETCDF_CFLAGS)
target_compile_options(benchmark PRIVATE ${NETCDF_CFLAGS})

# Add NetCDF linker flags
separate_arguments(NETCDF_FLIBS)
target_link_libraries(benchmark PRIVATE ${NETCDF_FLIBS} oasim)