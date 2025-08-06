cmake_minimum_required(VERSION 3.0)

add_subdirectory(f90)

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/f90)
    