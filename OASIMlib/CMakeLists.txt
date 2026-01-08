cmake_minimum_required(VERSION 3.0)

# ---------------------------------------------------------
# NVTX wrapper module
# ---------------------------------------------------------
# add_library(nvtx_wrapper nvtx_tools.f90)
# target_include_directories(nvtx_wrapper PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

# # Find CUDA toolkit (for libnvToolsExt)
# find_package(CUDAToolkit REQUIRED)
# target_link_libraries(nvtx_wrapper PUBLIC ${CUDAToolkit_NVTOOLSEXT_LIBRARIES})
# ---------------------------------------------------------
# Profiling flags
# ---------------------------------------------------------
# set(PROFILING_FLAGS "-pg")  
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# ---------------------------------------------------------
# Common object library
# ---------------------------------------------------------
add_library(common OBJECT common.f90)
set_property(TARGET common PROPERTY POSITION_INDEPENDENT_CODE ON)
target_compile_options(common PRIVATE ${PROFILING_FLAGS})
# ---------------------------------------------------------
# Main OASIM library
# ---------------------------------------------------------
add_library(oasim STATIC 
    oasim.f90 
    config.f90
    tables.f90
    matrix_read.f90
    # slingo_standalone.f90
    # slingo.f90
    # clrtrans.f90
    # clrtrans_standalone.f90
    # light_standalone.f90
    device_modules.f90
    # light.f90
    ocalbedo.f90
    sunmod.f90
    sfcirr.f90
    sfcsolz.f90
    monrad.f90
)

target_compile_options(oasim PRIVATE ${PROFILING_FLAGS})
target_link_libraries(oasim PUBLIC yaml common)
target_compile_options(oasim PRIVATE -cpp)
# target_link_libraries(oasim PRIVATE nvtx_wrapper)

# ---------------------------------------------------------
# Python interface library
# ---------------------------------------------------------
add_library(oasim-py SHARED oasim-py.f90)
target_compile_options(oasim-py PRIVATE ${PROFILING_FLAGS})
target_link_libraries(oasim-py PUBLIC oasim)

# if(USE_MPI)
#     target_link_libraries(oasim PUBLIC MPI::MPI_Fortran)
#     target_compile_definitions(oasim PUBLIC USE_MPI)
# endif()

# ---------------------------------------------------------
# Subdirectories
# ---------------------------------------------------------
add_subdirectory(fortran-yaml)

