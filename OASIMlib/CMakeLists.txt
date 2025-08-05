cmake_minimum_required(VERSION 3.0)

project(oasim LANGUAGES Fortran)


add_library(common OBJECT common.f90)
set_property(TARGET common PROPERTY POSITION_INDEPENDENT_CODE ON)
# NOTE: I changed the linking from SHARED TO STATIC to enable profiling 
add_library(${PROJECT_NAME} STATIC oasim.f90 
                         config.f90
                         tables.f90
                         matrix_read.f90
                         slingo.f90
                         clrtrans.f90
                         light.f90
                         ocalbedo.f90
                         sunmod.f90
                         sfcirr.f90
                         sfcsolz.f90
                         monrad.f90)
target_link_libraries(${PROJECT_NAME} PUBLIC yaml common)
# NOTE: I changed the linking from SHARED TO STATIC to enable profiling 
add_library(oasim-py STATIC oasim-py.f90)
target_link_libraries(oasim-py PUBLIC oasim)

add_subdirectory(fortran-yaml)
