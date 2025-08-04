cmake_minimum_required(VERSION 3.0)


add_executable(test test.f90)

set_target_properties(test PROPERTIES SUFFIX ".x")

separate_arguments(NETCDF_CFLAGS)
target_compile_options(test PRIVATE ${NETCDF_CFLAGS})

target_link_libraries(test PRIVATE ${NETCDF_FLIBS} oasim)