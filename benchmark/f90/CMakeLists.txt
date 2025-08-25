cmake_minimum_required(VERSION 3.0)

execute_process(COMMAND nf-config --fflags
    OUTPUT_VARIABLE NETCDF_FFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nf-config --flibs
    OUTPUT_VARIABLE NETCDF_FLIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

separate_arguments(NETCDF_FFLAGS NATIVE_COMMAND ${NETCDF_FFLAGS})
separate_arguments(NETCDF_FLIBS  NATIVE_COMMAND ${NETCDF_FLIBS})

if(CMAKE_Fortran_COMPILER MATCHES "ifx")
    message(STATUS "Intel Fortran compiler detected: appending imf")
    list(APPEND NETCDF_FLIBS imf)
endif()



add_library(NetCDF::Fortran INTERFACE IMPORTED)
target_compile_options(NetCDF::Fortran INTERFACE ${NETCDF_FFLAGS})
target_link_libraries(NetCDF::Fortran INTERFACE ${NETCDF_FLIBS})

add_executable(benchmark benchmark.f90)
set_target_properties(benchmark PROPERTIES SUFFIX ".x")

target_link_libraries(benchmark PRIVATE NetCDF::Fortran oasim)
