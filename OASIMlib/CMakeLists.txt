cmake_minimum_required(VERSION 3.0)
 
add_library(common OBJECT common.f90)
set_property(TARGET common PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(oasim SHARED oasim.f90 
                         config.f90
                         tables.f90
                         matrix_read.f90
                         slingo_standalone.f90
                        #  slingo.f90
                        #  clrtrans.f90
                        clrtrans_standalone.f90
                        light_standalone.f90
                        #  light.f90
                         ocalbedo.f90
                         sunmod.f90
                         sfcirr.f90
                         sfcsolz.f90
                         monrad.f90)
target_link_libraries(oasim PUBLIC yaml common)

add_library(oasim-py SHARED oasim-py.f90)
target_link_libraries(oasim-py PUBLIC oasim)

add_subdirectory(fortran-yaml)
