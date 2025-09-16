cmake_minimum_required(VERSION 3.0)

execute_process(COMMAND nf-config --fflags
    OUTPUT_VARIABLE NETCDF_FFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND nf-config --flibs
#     OUTPUT_VARIABLE NETCDF_FLIBS
#     OUTPUT_STRIP_TRAILING_WHITESPACE)
cmake_minimum_required(VERSION 3.0)

# Hard-coded NetCDF Fortran libraries (without -lm)
set(NETCDF_FLIBS
    -L/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/oneapi-2024.1.0/netcdf-fortran-4.6.1-gtz562rmr2nk57j7w2g32oniq4xokxak/lib
    -lnetcdff
    -lnetcdf
    -lnetcdf
    -limf
    # -lm  <-- removed to test if it fixes the linker
)
separate_arguments(NETCDF_FFLAGS NATIVE_COMMAND ${NETCDF_FFLAGS})
# separate_arguments(NETCDF_FLIBS  NATIVE_COMMAND ${NETCDF_FLIBS})


add_library(NetCDF::Fortran INTERFACE IMPORTED)
target_compile_options(NetCDF::Fortran INTERFACE ${NETCDF_FFLAGS})
target_link_libraries(NetCDF::Fortran INTERFACE ${NETCDF_FLIBS})

add_executable(benchmark benchmark.f90)
set_target_properties(benchmark PROPERTIES SUFFIX ".x")

target_link_libraries(benchmark PRIVATE NetCDF::Fortran oasim)