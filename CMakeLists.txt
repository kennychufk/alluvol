cmake_minimum_required (VERSION 3.20.0)
project(alluvol
  LANGUAGES CXX C
)

option(USE_SYSTEM_OPENVDB "Use system-installed OpenVDB instead of downloading a standalone OpenVDB" OFF)

include(FetchContent)

if (USE_SYSTEM_OPENVDB)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  find_package(OpenVDB REQUIRED)
else()
  FetchContent_Declare(
    openvdb
    GIT_REPOSITORY https://github.com/AcademySoftwareFoundation/openvdb.git
    GIT_TAG        v11.0.0
  )
  FetchContent_MakeAvailable(openvdb)
endif()

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.13.1
)
FetchContent_MakeAvailable(pybind11)

list(APPEND lib_alluvol_SOURCES
  src/io.cpp
  src/particle_list.cpp
  src/level_set.cpp
  src/geometry.cpp
)

pybind11_add_module(_ext
  src/pyalluvol.cpp
  ${lib_alluvol_SOURCES}
)

target_compile_definitions(_ext
  PRIVATE VERSION_INFO=${ALLUVOL_VERSION_INFO})

target_include_directories(_ext
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

target_link_libraries(_ext PRIVATE openvdb)
