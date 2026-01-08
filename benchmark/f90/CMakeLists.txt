cmake_minimum_required(VERSION 3.0)

option(ENABLE_PROFILING "Enable profiling with gprof" ON)

# ---------------------------------------------------------
# Define profiling flags
# ---------------------------------------------------------
if(ENABLE_PROFILING)
    # set(PROFILING_FLAGS "-pg")  
    message(STATUS "Enabling profiling flags: ${PROFILING_FLAGS}")
else()
    set(PROFILING_FLAGS "")
endif()

# ---------------------------------------------------------
# Find NetCDF Fortran
# ---------------------------------------------------------
find_package(NetCDF QUIET COMPONENTS Fortran)

if(NOT NetCDF_FOUND)
    message(STATUS "NetCDF package not found — using nf-config fallback")

    execute_process(COMMAND nf-config --fflags
        OUTPUT_VARIABLE NETCDF_FFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND nf-config --flibs
        OUTPUT_VARIABLE NETCDF_FLIBS
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    separate_arguments(NETCDF_FFLAGS)
    separate_arguments(NETCDF_FLIBS)

    add_library(NetCDF::Fortran INTERFACE IMPORTED)
    target_compile_options(NetCDF::Fortran INTERFACE ${NETCDF_FFLAGS})
    target_link_libraries(NetCDF::Fortran INTERFACE ${NETCDF_FLIBS})
endif()

# ---------------------------------------------------------
# Benchmark executable
# ---------------------------------------------------------
add_executable(benchmark benchmark.f90)
set_target_properties(benchmark PROPERTIES SUFFIX ".x")

# Apply profiling flags
target_compile_options(benchmark PRIVATE ${PROFILING_FLAGS})
target_link_libraries(benchmark PRIVATE NetCDF::Fortran oasim)
# target_link_options(benchmark PRIVATE -pg)
if(USE_MPI)
    target_link_libraries(benchmark PUBLIC MPI::MPI_Fortran)
    target_compile_definitions(benchmark PUBLIC USE_MPI)
endif()
