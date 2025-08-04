cmake_minimum_required(VERSION 3.0)

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