cmake_minimum_required(VERSION 3.20.0)
project(alluvion
  LANGUAGES CXX C
)
# set (CMAKE_CUDA_STANDARD 17)

find_package(OpenVDB REQUIRED)

# add_subdirectory(lib)
# add_subdirectory(tests)
add_subdirectory(src)
