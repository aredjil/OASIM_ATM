cmake_minimum_required(VERSION 3.0)

add_subdirectory(f90)
# Copy the data files to the build directory at configuration time 
# The test excutables require config.yaml and data directory to run
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY python/oasim.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python)

file(COPY config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/f90)
file(COPY config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python)
